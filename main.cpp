#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <oneapi/tbb/global_control.h>

using namespace std;

static const double G = 6.67e-11;   // Gravitational constant
static const double THETA = 0.5;    // Barnes-Hut opening angle
static const int MAX_DEPTH = 64;    // Stop subdividing (bodies at the same position)

struct Point {
    double x;
    double y;

    Point(double _x = 0, double _y = 0) : x(_x), y(_y) {}
};

struct Body {
    Point pos;
    double Fx, Fy;   // Net force accumulated during the current step
    double Vx, Vy, M;
    string S;

    Body(Point _pos, double _Vx, double _Vy, double _M, string _S)
        : pos(_pos), Fx(0), Fy(0), Vx(_Vx), Vy(_Vy), M(_M), S(_S) {}
};

class BHTree {
    Point topLeft;
    Point botRight;
    int depth;

    vector<Body*> bodies;      // Bodies stored in this node (only in leaves)
    BHTree* child[4];          // NW, NE, SW, SE

    double totalMass;
    Point centerOfMass;

    bool isLeaf() const {
        return !child[0] && !child[1] && !child[2] && !child[3];
    }

    int quadrant(const Point& p) const {
        double midX = (topLeft.x + botRight.x) / 2;
        double midY = (topLeft.y + botRight.y) / 2;
        bool west  = p.x <= midX;
        bool north = p.y <= midY;
        if (north) return west ? 0 : 1;
        return west ? 2 : 3;
    }

    void insertIntoChild(Body* body) {
        int q = quadrant(body->pos);
        if (!child[q]) {
            double midX = (topLeft.x + botRight.x) / 2;
            double midY = (topLeft.y + botRight.y) / 2;
            Point tl(q % 2 == 0 ? topLeft.x : midX, q < 2 ? topLeft.y : midY);
            Point br(q % 2 == 0 ? midX : botRight.x, q < 2 ? midY : botRight.y);
            child[q] = new BHTree(tl, br, depth + 1);
        }
        child[q]->insert(body);
    }

    static void addForce(Body* body, const Point& from, double mass) {
        double dx = from.x - body->pos.x;
        double dy = from.y - body->pos.y;
        double r = sqrt(dx * dx + dy * dy);
        if (r == 0) {
            return;
        }
        double F = (G * mass * body->M) / (r * r);
        body->Fx += F * (dx / r);
        body->Fy += F * (dy / r);
    }

public:
    BHTree(Point topL, Point botR, int _depth = 0)
        : topLeft(topL), botRight(botR), depth(_depth), totalMass(0) {
        child[0] = child[1] = child[2] = child[3] = nullptr;
    }

    ~BHTree() {
        for (BHTree* c : child) {
            delete c;
        }
    }

    // The tree owns raw pointers to its children, so it must not be copied
    BHTree(const BHTree&) = delete;
    BHTree& operator=(const BHTree&) = delete;

    // Check if the region of this node contains the point
    bool contains(const Point& p) const {
        return p.x >= topLeft.x && p.x <= botRight.x
            && p.y >= topLeft.y && p.y <= botRight.y;
    }

    // Insert a body. Bodies outside the universe are ignored.
    void insert(Body* body) {
        if (body == nullptr || !contains(body->pos)) {
            return;
        }

        if (isLeaf()) {
            if (bodies.empty() || depth >= MAX_DEPTH) {
                bodies.push_back(body);
                return;
            }
            // Leaf already has a body: split it and push the old body down
            vector<Body*> old;
            old.swap(bodies);
            for (Body* b : old) {
                insertIntoChild(b);
            }
        }
        insertIntoChild(body);
    }

    // Compute total mass and center of mass of every node (children in parallel)
    void calculateCenterOfMass() {
        totalMass = 0;
        double cx = 0, cy = 0;

        if (isLeaf()) {
            for (Body* b : bodies) {
                totalMass += b->M;
                cx += b->M * b->pos.x;
                cy += b->M * b->pos.y;
            }
        } else {
            tbb::parallel_invoke(
                [this]() { if (child[0]) child[0]->calculateCenterOfMass(); },
                [this]() { if (child[1]) child[1]->calculateCenterOfMass(); },
                [this]() { if (child[2]) child[2]->calculateCenterOfMass(); },
                [this]() { if (child[3]) child[3]->calculateCenterOfMass(); }
            );
            for (BHTree* c : child) {
                if (c) {
                    totalMass += c->totalMass;
                    cx += c->totalMass * c->centerOfMass.x;
                    cy += c->totalMass * c->centerOfMass.y;
                }
            }
        }

        if (totalMass > 0) {
            centerOfMass = Point(cx / totalMass, cy / totalMass);
        } else {
            centerOfMass = Point((topLeft.x + botRight.x) / 2, (topLeft.y + botRight.y) / 2);
        }
    }

    // Accumulate the force that this node exerts on the body
    void netForce(Body* body) const {
        if (totalMass == 0) {
            return;
        }

        if (isLeaf()) {
            for (Body* other : bodies) {
                if (other != body) {
                    addForce(body, other->pos, other->M);
                }
            }
            return;
        }

        double dx = centerOfMass.x - body->pos.x;
        double dy = centerOfMass.y - body->pos.y;
        double r = sqrt(dx * dx + dy * dy);
        double width = botRight.x - topLeft.x;

        // Far enough away (and not containing the body): treat the node as one mass
        if (!contains(body->pos) && r > 0 && width / r < THETA) {
            addForce(body, centerOfMass, totalMass);
        } else {
            for (BHTree* c : child) {
                if (c) c->netForce(body);
            }
        }
    }
};

// Update velocity and position from the accumulated force, then reset it
static void newPosition(Body& body, double dt) {
    body.Vx += dt * body.Fx / body.M;
    body.Vy += dt * body.Fy / body.M;

    body.pos.x += dt * body.Vx;
    body.pos.y += dt * body.Vy;

    body.Fx = 0;
    body.Fy = 0;
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 5) {
        cout << "Usage: " << argv[0] << " input_file.txt iterations [threads] [dt]" << endl;
        return 1;
    }

    string inputFileName = argv[1];
    int numIterations = atoi(argv[2]);
    size_t numThreads = argc > 3 ? (size_t)atoi(argv[3]) : max(1u, thread::hardware_concurrency());
    double dt = argc > 4 ? atof(argv[4]) : 1.0;

    if (numIterations < 0 || numThreads < 1 || dt <= 0) {
        cerr << "Error: iterations must be >= 0, threads >= 1 and dt > 0." << endl;
        return 1;
    }

    oneapi::tbb::global_control gc(oneapi::tbb::global_control::max_allowed_parallelism, numThreads);

    // Open input file
    ifstream inputFile(inputFileName);
    if (!inputFile) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }

    // Read simulation parameters from the input file
    int N;      // Number of bodies
    double R;   // Radius of the universe (positions are in [-R, R])
    if (!(inputFile >> N >> R)) {
        cerr << "Error: Invalid input file header." << endl;
        return 1;
    }

    // Read body data from the input file
    vector<Body> bodies;
    bodies.reserve(N);
    for (int i = 0; i < N; ++i) {
        double X, Y, Vx, Vy, M;
        string S;
        if (!(inputFile >> X >> Y >> Vx >> Vy >> M >> S)) {
            cerr << "Error: Expected " << N << " bodies but could only read " << i << "." << endl;
            return 1;
        }
        bodies.push_back(Body(Point(X, Y), Vx, Vy, M, S));
    }
    inputFile.close();

    // Perform simulation iterations
    for (int iter = 0; iter < numIterations; ++iter) {
        // Build a new tree for the current positions. Insertion is sequential
        // because insert() modifies shared nodes and is not thread-safe.
        BHTree tree(Point(-R, -R), Point(R, R));
        for (Body& body : bodies) {
            tree.insert(&body);
        }

        tree.calculateCenterOfMass();

        // Calculate forces in parallel (each task only writes its own body)
        tbb::parallel_for(tbb::blocked_range<size_t>(0, bodies.size()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                tree.netForce(&bodies[i]);
            }
        });

        // Update positions and velocities in parallel
        tbb::parallel_for(tbb::blocked_range<size_t>(0, bodies.size()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                newPosition(bodies[i], dt);
            }
        });
    }

    // Write output data in the same format as the input
    ofstream outputFile("output.txt");
    if (!outputFile) {
        cerr << "Error: Could not open output file." << endl;
        return 1;
    }

    outputFile << N << endl;
    outputFile << scientific << setprecision(4) << R << endl;
    for (const Body& body : bodies) {
        outputFile << body.pos.x << " " << body.pos.y << " " << body.Vx << " " << body.Vy << " "
                   << body.M << " " << body.S << endl;
    }
    outputFile.close();

    return 0;
}
