import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.ForkJoinPool;



class Point {
    public int x;
    public int y;

    public Point(int _x, int _y) {
        x = _x;
        y = _y;
    }

    public Point() {
        x = 0;
        y = 0;
    }
}

class Body {
    public Point pos;
    public double X, Y, Vx, Vy, M;
    public String S;

    public Body(Point _pos, double _X, double _Y, double _Vx, double _Vy, double _M, String _S) {
        pos = _pos;
        X = _X;
        Y = _Y;
        Vx = _Vx;
        Vy = _Vy;
        M = _M;
        S = _S;
    }

    public Body() {
        X = 0;
        Y = 0;
        Vx = 0;
        Vy = 0;
        M = 0;
        S = "";
    }
}

class BHTree {
    private Point topLeft;
    private Point botRight;

    private Body b;

    private BHTree NW;
    private BHTree NE;
    private BHTree SW;
    private BHTree SE;

    private double totalMass;
    private Point centerOfMass;

    public BHTree() {
        topLeft = new Point(0, 0);
        botRight = new Point(0, 0);
        b = null;
        NW = null;
        NE = null;
        SW = null;
        SE = null;
    }

    public BHTree(Point topL, Point botR) {
        topLeft = topL;
        botRight = botR;
        b = null;
        NW = null;
        NE = null;
        SW = null;
        SE = null;
    }

    public Point getCenterOfMass() {
        return centerOfMass;
    }

    public void insert(Body body) {
        if (body == null) {
            return;
        }

        if (!length(body.pos)) {
            return;
        }

        if (Math.abs(topLeft.x - botRight.x) <= 1 && Math.abs(topLeft.y - botRight.y) <= 1) {
            if (b == null) {
                b = body;
            }
            return;
        }

        if ((topLeft.x + botRight.x) / 2 >= body.pos.x) {
            if ((topLeft.y + botRight.y) / 2 >= body.pos.y) {
                if (NW == null)
                    NW = new BHTree(
                            new Point(topLeft.x, topLeft.y),
                            new Point((topLeft.x + botRight.x) / 2,
                                    (topLeft.y + botRight.y) / 2));
                NW.insert(body);
            } else {
                if (SW == null)
                    SW = new BHTree(
                            new Point(topLeft.x,
                                    (topLeft.y + botRight.y) / 2),
                            new Point((topLeft.x + botRight.x) / 2,
                                    botRight.y));
                SW.insert(body);
            }
        } else {
            if ((topLeft.y + botRight.y) / 2 >= body.pos.y) {
                if (NE == null)
                    NE = new BHTree(
                            new Point((topLeft.x + botRight.x) / 2,
                                    topLeft.y),
                            new Point(botRight.x,
                                    (topLeft.y + botRight.y) / 2));
                NE.insert(body);
            } else {
                if (SE == null)
                    SE = new BHTree(
                            new Point((topLeft.x + botRight.x) / 2,
                                    (topLeft.y + botRight.y) / 2),
                            new Point(botRight.x, botRight.y));
                SE.insert(body);
            }
        }
    }

    private boolean length(Point p) {
        return (p.x >= topLeft.x && p.x <= botRight.x
                && p.y >= topLeft.y && p.y <= botRight.y);
    }

    public void calculateCenterOfMassParallel() {
        new CalculateCenterOfMassTask(this).invoke();
    }

    private static class CalculateCenterOfMassTask extends RecursiveAction {
        private static final int THRESHOLD = 100; // adjust as needed
        private BHTree node;

        public CalculateCenterOfMassTask(BHTree node) {
            this.node = node;
        }

        @Override
        protected void compute() {
            if (node.NW == null && node.NE == null && node.SW == null && node.SE == null) {
                if (node.b != null) {
                    node.totalMass = node.b.M;
                    node.centerOfMass = node.b.pos;
                } else {
                    node.totalMass = 0;
                    node.centerOfMass = new Point((node.topLeft.x + node.botRight.x) / 2, (node.topLeft.y + node.botRight.y) / 2);
                }
                return;
            }

            if (node.NW != null) node.NW.calculateCenterOfMassParallel();
            if (node.NE != null) node.NE.calculateCenterOfMassParallel();
            if (node.SW != null) node.SW.calculateCenterOfMassParallel();
            if (node.SE != null) node.SE.calculateCenterOfMassParallel();

            node.totalMass = 0;
            node.centerOfMass = new Point(0, 0);
            if (node.NW != null) {
                node.totalMass += node.NW.totalMass;
                node.centerOfMass.x += node.NW.totalMass * node.NW.centerOfMass.x;
                node.centerOfMass.y += node.NW.totalMass * node.NW.centerOfMass.y;
            }
            if (node.NE != null) {
                node.totalMass += node.NE.totalMass;
                node.centerOfMass.x += node.NE.totalMass * node.NE.centerOfMass.x;
                node.centerOfMass.y += node.NE.totalMass * node.NE.centerOfMass.y;
            }
            if (node.SW != null) {
                node.totalMass += node.SW.totalMass;
                node.centerOfMass.x += node.SW.totalMass * node.SW.centerOfMass.x;
                node.centerOfMass.y += node.SW.totalMass * node.SW.centerOfMass.y;
            }
            if (node.SE != null) {
                node.totalMass += node.SE.totalMass;
                node.centerOfMass.x += node.SE.totalMass * node.SE.centerOfMass.x;
                node.centerOfMass.y += node.SE.totalMass * node.SE.centerOfMass.y;
            }

            if (node.totalMass > 0) {
                node.centerOfMass.x /= node.totalMass;
                node.centerOfMass.y /= node.totalMass;
            }
        }
    }

    public void netForce(Body body) {
        if (!length(body.pos) || (Math.abs(topLeft.x - botRight.x) <= 1 && Math.abs(topLeft.y - botRight.y) <= 1) || b == body) {
            return;
        }

        double dx = centerOfMass.x - body.pos.x;
        double dy = centerOfMass.y - body.pos.y;
        double r = Math.sqrt(dx * dx + dy * dy);

        if (r > Math.min(Math.abs(topLeft.x - botRight.x), Math.abs(topLeft.y - botRight.y))) {
            double F = (6.67e-11 * totalMass * body.M) / (r * r);
            double Fx = F * (dx / r);
            double Fy = F * (dy / r);
            body.X += Fx;
            body.Y += Fy;
        } else {
            if (NW != null) NW.netForce(body);
            if (NE != null) NE.netForce(body);
            if (SW != null) SW.netForce(body);
            if (SE != null) SE.netForce(body);
        }
    }

    public void newPosition(Body body) {
        body.Vx += body.X / body.M;
        body.Vy += body.Y / body.M;
        body.pos.x += body.Vx;
        body.pos.y += body.Vy;
        body.X = 0;
        body.Y = 0;
    }

    public void computeNetForcesParallel(List<Body> bodies) {
        bodies.parallelStream().forEach(this::netForce);
    }
}

public class BarnesHutSimulation {
    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Usage: BarnesHutSimulation input_file.txt Iter");
            return;
        }

        String inputFileName = args[0];
        int numIterations = Integer.parseInt(args[1]);

        int numThreads = 4; // Adjust as needed
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", String.valueOf(numThreads));

        List<Body> bodies = new ArrayList<>();
        double R = 0.0; // Default value

        try (Scanner scanner = new Scanner(new File(inputFileName))) {
            int N = scanner.nextInt();
            R = scanner.nextDouble(); // Read simulation parameter

            for (int i = 0; i < N; i++) {
                double X = scanner.nextDouble();
                double Y = scanner.nextDouble();
                double Vx = scanner.nextDouble();
                double Vy = scanner.nextDouble();
                double M = scanner.nextDouble();
                String S = scanner.next();
                bodies.add(new Body(new Point((int) X, (int) Y), X, Y, Vx, Vy, M, S));
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error: Could not open input file.");
            System.exit(1);
        }

        ForkJoinPool pool = new ForkJoinPool(); // Create a ForkJoinPool

        // Define a task for parallel execution
        RecursiveAction simulationTask = new RecursiveAction() {
            @Override
            protected void compute() {
                for (int iter = 0; iter < numIterations; ++iter) {
                    BHTree bhtree = new BHTree(new Point(0, 0), new Point(100, 100)); // Create a new instance
                    // Construct the Barnes-Hut tree
                    for (Body body : bodies) {
                        bhtree.insert(body);
                    }
                    bhtree.calculateCenterOfMassParallel(); // Use the parallel version

                    // Calculate forces using Barnes-Hut algorithm in parallel
                    bhtree.computeNetForcesParallel(bodies);

                    // Update positions and velocities
                    for (Body body : bodies) {
                        bhtree.newPosition(body);
                    }
                }
            }
        };

        // Execute the simulation task in parallel
        pool.invoke(simulationTask);

        // Write the simulation results to a file
        String outputFileName = "output.txt";
        try (FileWriter writer = new FileWriter(outputFileName)) {
            writer.write(bodies.size() + "\n"); // Number of bodies

            // Write simulation parameter R
            writer.write(R + "\n");

            // Write body information
            for (Body body : bodies) {
                writer.write(body.pos.x + " " + body.pos.y + " " + body.Vx + " " + body.Vy + " " + body.M + " " + body.S + "\n");
            }
        } catch (IOException e) {
            System.err.println("Error: Could not open or write to output file.");
            System.exit(1);
        }
    }
}