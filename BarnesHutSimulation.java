import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

class Point {
    public double x;
    public double y;

    public Point(double _x, double _y) {
        x = _x;
        y = _y;
    }

    public Point() {
        this(0, 0);
    }
}

class Body {
    public Point pos;
    public double Fx, Fy; // Net force accumulated during the current step
    public double Vx, Vy, M;
    public String S;

    public Body(Point _pos, double _Vx, double _Vy, double _M, String _S) {
        pos = _pos;
        Vx = _Vx;
        Vy = _Vy;
        M = _M;
        S = _S;
    }
}

class BHTree {
    private static final double G = 6.67e-11;  // Gravitational constant
    private static final double THETA = 0.5;   // Barnes-Hut opening angle
    private static final int MAX_DEPTH = 64;   // Stop subdividing (bodies at the same position)

    private final Point topLeft;
    private final Point botRight;
    private final int depth;

    private final List<Body> bodies = new ArrayList<>(); // Bodies stored in this node (only in leaves)
    private final BHTree[] child = new BHTree[4];         // NW, NE, SW, SE

    private double totalMass;
    private Point centerOfMass = new Point();

    public BHTree(Point topL, Point botR) {
        this(topL, botR, 0);
    }

    private BHTree(Point topL, Point botR, int depth) {
        topLeft = topL;
        botRight = botR;
        this.depth = depth;
    }

    private boolean isLeaf() {
        return child[0] == null && child[1] == null && child[2] == null && child[3] == null;
    }

    // Check if the region of this node contains the point
    public boolean contains(Point p) {
        return p.x >= topLeft.x && p.x <= botRight.x
                && p.y >= topLeft.y && p.y <= botRight.y;
    }

    private int quadrant(Point p) {
        double midX = (topLeft.x + botRight.x) / 2;
        double midY = (topLeft.y + botRight.y) / 2;
        boolean west = p.x <= midX;
        boolean north = p.y <= midY;
        if (north) return west ? 0 : 1;
        return west ? 2 : 3;
    }

    private void insertIntoChild(Body body) {
        int q = quadrant(body.pos);
        if (child[q] == null) {
            double midX = (topLeft.x + botRight.x) / 2;
            double midY = (topLeft.y + botRight.y) / 2;
            Point tl = new Point(q % 2 == 0 ? topLeft.x : midX, q < 2 ? topLeft.y : midY);
            Point br = new Point(q % 2 == 0 ? midX : botRight.x, q < 2 ? midY : botRight.y);
            child[q] = new BHTree(tl, br, depth + 1);
        }
        child[q].insert(body);
    }

    // Insert a body. Bodies outside the universe are ignored.
    public void insert(Body body) {
        if (body == null || !contains(body.pos)) {
            return;
        }

        if (isLeaf()) {
            if (bodies.isEmpty() || depth >= MAX_DEPTH) {
                bodies.add(body);
                return;
            }
            // Leaf already has a body: split it and push the old body down
            List<Body> old = new ArrayList<>(bodies);
            bodies.clear();
            for (Body b : old) {
                insertIntoChild(b);
            }
        }
        insertIntoChild(body);
    }

    public void calculateCenterOfMassParallel() {
        new CalculateCenterOfMassTask(this).invoke();
    }

    // Computes total mass / center of mass; the four children are processed in parallel
    private static class CalculateCenterOfMassTask extends RecursiveAction {
        private final BHTree node;

        CalculateCenterOfMassTask(BHTree node) {
            this.node = node;
        }

        @Override
        protected void compute() {
            double total = 0, cx = 0, cy = 0;

            if (node.isLeaf()) {
                for (Body b : node.bodies) {
                    total += b.M;
                    cx += b.M * b.pos.x;
                    cy += b.M * b.pos.y;
                }
            } else {
                List<CalculateCenterOfMassTask> tasks = new ArrayList<>();
                for (BHTree c : node.child) {
                    if (c != null) tasks.add(new CalculateCenterOfMassTask(c));
                }
                invokeAll(tasks);

                for (BHTree c : node.child) {
                    if (c != null) {
                        total += c.totalMass;
                        cx += c.totalMass * c.centerOfMass.x;
                        cy += c.totalMass * c.centerOfMass.y;
                    }
                }
            }

            node.totalMass = total;
            if (total > 0) {
                node.centerOfMass = new Point(cx / total, cy / total);
            } else {
                node.centerOfMass = new Point((node.topLeft.x + node.botRight.x) / 2,
                        (node.topLeft.y + node.botRight.y) / 2);
            }
        }
    }

    private static void addForce(Body body, Point from, double mass) {
        double dx = from.x - body.pos.x;
        double dy = from.y - body.pos.y;
        double r = Math.sqrt(dx * dx + dy * dy);
        if (r == 0) {
            return;
        }
        double F = (G * mass * body.M) / (r * r);
        body.Fx += F * (dx / r);
        body.Fy += F * (dy / r);
    }

    // Accumulate the force that this node exerts on the body
    public void netForce(Body body) {
        if (totalMass == 0) {
            return;
        }

        if (isLeaf()) {
            for (Body other : bodies) {
                if (other != body) {
                    addForce(body, other.pos, other.M);
                }
            }
            return;
        }

        double dx = centerOfMass.x - body.pos.x;
        double dy = centerOfMass.y - body.pos.y;
        double r = Math.sqrt(dx * dx + dy * dy);
        double width = botRight.x - topLeft.x;

        // Far enough away (and not containing the body): treat the node as one mass
        if (!contains(body.pos) && r > 0 && width / r < THETA) {
            addForce(body, centerOfMass, totalMass);
        } else {
            for (BHTree c : child) {
                if (c != null) c.netForce(body);
            }
        }
    }

    // Update velocity and position from the accumulated force, then reset it
    public static void newPosition(Body body, double dt) {
        body.Vx += dt * body.Fx / body.M;
        body.Vy += dt * body.Fy / body.M;
        body.pos.x += dt * body.Vx;
        body.pos.y += dt * body.Vy;
        body.Fx = 0;
        body.Fy = 0;
    }
}

public class BarnesHutSimulation {
    public static void main(String[] args) {
        if (args.length < 2 || args.length > 4) {
            System.out.println("Usage: java BarnesHutSimulation input_file.txt iterations [threads] [dt]");
            return;
        }

        String inputFileName = args[0];
        int numIterations;
        int numThreads;
        double dt;
        try {
            numIterations = Integer.parseInt(args[1]);
            numThreads = args.length > 2 ? Integer.parseInt(args[2]) : Runtime.getRuntime().availableProcessors();
            dt = args.length > 3 ? Double.parseDouble(args[3]) : 1.0;
        } catch (NumberFormatException e) {
            System.err.println("Error: iterations and threads must be integers, dt a number.");
            return;
        }
        if (numIterations < 0 || numThreads < 1 || dt <= 0) {
            System.err.println("Error: iterations must be >= 0, threads >= 1 and dt > 0.");
            return;
        }

        List<Body> bodies = new ArrayList<>();
        int N;
        double R;

        try (Scanner scanner = new Scanner(new File(inputFileName))) {
            scanner.useLocale(Locale.US);
            N = scanner.nextInt();
            R = scanner.nextDouble();

            for (int i = 0; i < N; i++) {
                double X = scanner.nextDouble();
                double Y = scanner.nextDouble();
                double Vx = scanner.nextDouble();
                double Vy = scanner.nextDouble();
                double M = scanner.nextDouble();
                String S = scanner.next();
                bodies.add(new Body(new Point(X, Y), Vx, Vy, M, S));
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error: Could not open input file.");
            System.exit(1);
            return;
        } catch (NoSuchElementException e) { // also covers InputMismatchException
            System.err.println("Error: Invalid input file format.");
            System.exit(1);
            return;
        }

        final double universe = R;
        final double step = dt;

        // All parallel work (fork/join tasks and parallel streams) runs inside this pool,
        // so the number of threads is controlled by numThreads.
        ForkJoinPool pool = new ForkJoinPool(numThreads);
        try {
            pool.submit(() -> {
                for (int iter = 0; iter < numIterations; ++iter) {
                    // Build a new tree for the current positions (insertion is sequential)
                    BHTree bhtree = new BHTree(new Point(-universe, -universe), new Point(universe, universe));
                    for (Body body : bodies) {
                        bhtree.insert(body);
                    }

                    bhtree.calculateCenterOfMassParallel();

                    // Calculate forces in parallel (each task only writes its own body)
                    bodies.parallelStream().forEach(bhtree::netForce);

                    // Update positions and velocities in parallel
                    bodies.parallelStream().forEach(b -> BHTree.newPosition(b, step));
                }
            }).get();
        } catch (InterruptedException | ExecutionException e) {
            System.err.println("Error: Simulation failed: " + e.getMessage());
            System.exit(1);
        } finally {
            pool.shutdown();
        }

        // Write the simulation results in the same format as the input
        try (FileWriter writer = new FileWriter("output.txt")) {
            writer.write(N + "\n");
            writer.write(String.format(Locale.US, "%.4e%n", R));
            for (Body body : bodies) {
                writer.write(String.format(Locale.US, "%.4e %.4e %.4e %.4e %.4e %s%n",
                        body.pos.x, body.pos.y, body.Vx, body.Vy, body.M, body.S));
            }
        } catch (IOException e) {
            System.err.println("Error: Could not open or write to output file.");
            System.exit(1);
        }
    }
}
