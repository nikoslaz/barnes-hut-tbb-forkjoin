#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <sstream>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <oneapi/tbb/global_control.h>


using namespace tbb;
using namespace std;

struct Point {
public:
    int x;
    int y;

    Point(double _x, double _y) {
        x = _x;
        y = _y;
    }

    Point() {
        x = 0;
        y = 0;
    }
};

struct Body{
public:
    Point pos;
    double X, Y, Vx, Vy, M;
    string S;

    Body(Point _pos, double _X, double _Y, double _Vx, double _Vy, double _M, string _S){
        pos = _pos;
        X   = _X;
        Y   = _Y;
        Vx  = _Vx;
        Vy  = _Vy;
        M   = _M;
        S   = _S;
    }

    Body(){
        X   = 0;
        Y   = 0;
        Vx  = 0;
        Vy  = 0;
        M   = 0;
        S   = "";
    }

};


class BHTree {
    Point topLeft;
    Point botRight;

    Body* b;

    BHTree* NW;
    BHTree* NE;
    BHTree* SW;
    BHTree* SE;

    double totalMass;
    Point centerOfMass;

public:

    BHTree() {
        topLeft = Point(0, 0);
        botRight = Point(0, 0);
        b = NULL;
        NW = NULL;
        NE = NULL;
        SW = NULL;
        SE = NULL;
    }

    BHTree(Point topL, Point botR) {
        topLeft  = topL;
        botRight = botR;
        b = NULL;
        NW = NULL;
        NE = NULL;
        SW = NULL;
        SE = NULL;
    }

    Point getCenterOfMass() const {
        return centerOfMass;
    }

    void insert(Body*);
    bool length(Point);
    void calculateCenterOfMass();
    void netForce(Body*);
    void newPosition(Body*);
    void computeNetForcesParallel(vector<Body*>&);

};

// Insert node in the BHTreetree
void BHTree::insert(Body* body) {

    if(body == NULL){
        return;
    }

    if(!length(body->pos)){
        return;
    }

    if(abs(topLeft.x - botRight.x) <=1
       && abs(topLeft.y - botRight.y) <=1) {
            if(b == NULL){
                b = body;
            }
            return;
    }

    if ((topLeft.x + botRight.x) / 2 >= body->pos.x) {
        // Indicates topLeftTree
        if ((topLeft.y + botRight.y) / 2 >= body->pos.y) {
            if (NW == NULL)
                NW = new BHTree(
                    Point(topLeft.x, topLeft.y),
                    Point((topLeft.x + botRight.x) / 2,
                          (topLeft.y + botRight.y) / 2));
            NW->insert(body);
        }
 
        // Indicates botLeftTree
        else {
            if (SW == NULL)
                SW = new BHTree(
                    Point(topLeft.x,
                          (topLeft.y + botRight.y) / 2),
                    Point((topLeft.x + botRight.x) / 2,
                          botRight.y));
            SW->insert(body);
        }
    }
    else {
        // Indicates topRightTree
        if ((topLeft.y + botRight.y) / 2 >= body->pos.y) {
            if (NE == NULL)
                NE = new BHTree(
                    Point((topLeft.x + botRight.x) / 2,
                          topLeft.y),
                    Point(botRight.x,
                          (topLeft.y + botRight.y) / 2));
            NE->insert(body);
        }
 
        // Indicates botRightTree
        else {
            if (SE == NULL)
                SE = new BHTree(
                    Point((topLeft.x + botRight.x) / 2,
                          (topLeft.y + botRight.y) / 2),
                    Point(botRight.x, botRight.y));
            SE->insert(body);
        }
    }

}

// Check if current BHTreetree contains the point
bool BHTree::length(Point p)
{
    return (p.x >= topLeft.x && p.x <= botRight.x
            && p.y >= topLeft.y && p.y <= botRight.y);
}

void BHTree::calculateCenterOfMass() {
    // if its a leaf
    if (NW == NULL && NE == NULL && SW == NULL && SE == NULL) {
        if (b != NULL) {
            totalMass = b->M;
            centerOfMass = b->pos;
        } else {
            totalMass = 0;
            centerOfMass = Point((topLeft.x + botRight.x) / 2, (topLeft.y + botRight.y) / 2);
        }
        return;
    }

    parallel_invoke(
        [this]() { if (NW != NULL) NW->calculateCenterOfMass(); },
        [this]() { if (NE != NULL) NE->calculateCenterOfMass(); },
        [this]() { if (SW != NULL) SW->calculateCenterOfMass(); },
        [this]() { if (SE != NULL) SE->calculateCenterOfMass(); }
    );

    // calculate total mass and the center of mass
    totalMass = 0;
    centerOfMass = Point(0, 0);
    if (NW != NULL) {
        totalMass += NW->totalMass;
        centerOfMass.x += NW->totalMass * NW->centerOfMass.x;
        centerOfMass.y += NW->totalMass * NW->centerOfMass.y;
    }
    if (NE != NULL) {
        totalMass += NE->totalMass;
        centerOfMass.x += NE->totalMass * NE->centerOfMass.x;
        centerOfMass.y += NE->totalMass * NE->centerOfMass.y;
    }
    if (SW != NULL) {
        totalMass += SW->totalMass;
        centerOfMass.x += SW->totalMass * SW->centerOfMass.x;
        centerOfMass.y += SW->totalMass * SW->centerOfMass.y;
    }
    if (SE != NULL) {
        totalMass += SE->totalMass;
        centerOfMass.x += SE->totalMass * SE->centerOfMass.x;
        centerOfMass.y += SE->totalMass * SE->centerOfMass.y;
    }

    if (totalMass > 0) {
        centerOfMass.x /= totalMass;
        centerOfMass.y /= totalMass;
    }
}

void BHTree::netForce(Body* body) {
    // Check if the current node can't subdivide further or contains the body
    if (!length(body->pos) || (abs(topLeft.x - botRight.x) <= 1 && abs(topLeft.y - botRight.y) <= 1) || b == body) {
        return;
    }

    // Calculate the distance between the center of mass of this node and the body
    double dx = centerOfMass.x - body->pos.x;
    double dy = centerOfMass.y - body->pos.y;
    double r = sqrt(dx * dx + dy * dy);

    // If the distance is greater than the size of the smaller region containing the body,
    // and the body is not contained in this node, calculate the force exerted on the body
    if (r > min(abs(topLeft.x - botRight.x), abs(topLeft.y - botRight.y))) {
        // Calculate the gravitational force using Newton's law of universal gravitation
        double F = (6.67e-11 * totalMass * body->M) / (r * r);

        // Calculate the force components in x and y directions
        double Fx = F * (dx / r);
        double Fy = F * (dy / r);

        // Update the net force exerted on the body
        body->X += Fx;
        body->Y += Fy;
    }
    // If the body is contained within this node or the distance is smaller than the size of the body,
    // recursively calculate the forces exerted by the children nodes
    else {
        if (NW != NULL) NW->netForce(body);
        if (NE != NULL) NE->netForce(body);
        if (SW != NULL) SW->netForce(body);
        if (SE != NULL) SE->netForce(body);
    }
}

void BHTree::newPosition(Body* body) {
    // Update velocity using the calculated net force
    body->Vx += body->X / body->M;
    body->Vy += body->Y / body->M;

    // Update position using the updated velocity
    body->pos.x += body->Vx;
    body->pos.y += body->Vy;

    // Reset net force for the next iteration
    body->X = 0;
    body->Y = 0;
}

void BHTree::computeNetForcesParallel(vector<Body*>& bodies) {
    parallel_for_each(bodies.begin(), bodies.end(), [&](Body* body) {
        netForce(body);
    });
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " input_file.txt Iter" << endl;
        return 1;
    }
    
    size_t numThreads = 2;
    oneapi::tbb::global_control gc(oneapi::tbb::global_control::max_allowed_parallelism, numThreads);

    // Read input file name from command line arguments
    string inputFileName = argv[1];
    int numIterations = atoi(argv[2]);

    // Open input file
    ifstream inputFile(inputFileName);
    if (!inputFile) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }

    // Read simulation parameters from the input file
    int N;      // Number of bodies
    double R;   // Size of the universe
    inputFile >> N >> R;

    // Read body data from the input file and store in a vector
    vector<Body> bodies;
    for (int i = 0; i < N; ++i) {
        double X, Y, Vx, Vy, M;
        string S;
        inputFile >> X >> Y >> Vx >> Vy >> M >> S;
        bodies.push_back(Body(Point(X, Y), X, Y, Vx, Vy, M, S));
    }

    // Close input file
    inputFile.close();

    // Use the same bounds as the tree that is rebuilt every iteration
    BHTree bhtree(Point(0, 0), Point(100, 100));

    // Create a vector of pointers to Body objects
    vector<Body*> bodyPointers;
    for (auto& body : bodies) {
        bodyPointers.push_back(&body);
    }

    // Perform simulation iterations
    for (int iter = 0; iter < numIterations; ++iter) {
        // Construct the Barnes-Hut tree. Insertion is sequential because
        // insert() modifies shared nodes and is not thread-safe.
        for (auto* body : bodyPointers) {
            bhtree.insert(body);
        }

        // Compute total mass / center of mass of every node (parallel_invoke inside)
        bhtree.calculateCenterOfMass();

        // Calculate forces using Barnes-Hut algorithm in parallel
        bhtree.computeNetForcesParallel(bodyPointers);

        // Update positions and velocities in parallel
        tbb::parallel_for(tbb::blocked_range<size_t>(0, bodies.size()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                bhtree.newPosition(&bodies[i]);
            }
        });

        // Clear the Barnes-Hut tree for the next iteration
        // If the tree is not cleared, the new bodies' positions
        // won't be correctly inserted into the tree
        bhtree = BHTree(Point(0, 0), Point(100, 100));
    }


    // Open output file
    ofstream outputFile("output.txt");
    if (!outputFile) {
        cerr << "Error: Could not open output file." << endl;
        return 1;
    }

    // Write output data to the output file
    outputFile << N << endl;
    outputFile << R << endl;
    for (const auto& body : bodies) {
        outputFile << body.pos.x << " " << body.pos.y << " " << body.Vx << " " << body.Vy << " " << body.M << " " << body.S << endl;
    }

    // Close output file
    outputFile.close();

    return 0;
}

