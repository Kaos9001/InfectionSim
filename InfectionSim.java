import java.awt.Color;

public class InfectionSim {
    private static final double HZ = 5;    // number of redraw events per clock tick
    private static final double HZ_SPREAD_CHECK = 5; // number of infection checks per clock tick
    private static double INFECTIVITY = 0.01; // Infectivity at max distance
    private static double BASE_INFECTIVITY = 1; // Global infectivity multiplier

    private MinPQ<Event> pq;          // the priority queue
    private double t  = 0.0;          // simulation clock time
    private Particle[] particles;     // the array of particles

    public Draw canvas;
    public Draw graph;

    private SIRCounter SIR_atual;


    /**
     * Initializes a system with the specified collection of particles.
     * The individual particles will be mutated during the simulation.
     *
     * @param  particles the array of particles
     */
    public InfectionSim(Particle[] particles) {
        this.particles = particles.clone();   // defensive copy
        this.canvas = new Draw("Simulacao");
        this.graph = new Draw("Grafico");
        this.graph.setYscale(0, 1);
        this.graph.setPenRadius(0.005);
        this.graph.setPenColor(Draw.RED);
        this.graph.setXscale(0, 20);
        this.graph.enableDoubleBuffering();
    }

    private class SIRCounter
    {
        public int S, I, R;

        public SIRCounter(int S, int I, int R)
        {
            this.S = S;
            this.I = I;
            this.R = R;
        }
    }


    // Adds future infection event after new particle is infected
    private void predict_zone_entry(Particle S, Particle I, double limit)
    {
        if (!S.inside_spread_zone)
        {
            double[] dts = S.timeToHit(I);

            if (dts[1] <= S.recovery_time && t + dts[1] <= limit)
            {
                pq.insert(new EnterSpreadZoneEvent(t + dts[1], S, I));
            }
        }
    }

    // updates priority queue with all new collision events for particle a
    private void predict(Particle a, double limit) {
        if (a == null) return;

        // particle-particle collisions
        for (int i = 0; i < particles.length; i++) 
        {
            double[] dts = a.timeToHit(particles[i]);

            if (t + dts[0] <= limit)
            {
                pq.insert(new CollisionEvent(t + dts[0], a, particles[i]));
            }

            if (t + dts[1] <= limit)
            {
                if (!particles[i].inside_spread_zone && a.infected && particles[i].susceptible())
                {
                    pq.insert(new EnterSpreadZoneEvent(t + dts[1], particles[i], a));
                }
                else if (!a.inside_spread_zone&& a.susceptible() && particles[i].infected)
                {
                    pq.insert(new EnterSpreadZoneEvent(t + dts[1], a, particles[i]));
                }
            }
        }

        // particle-wall collisions
        double dtX = a.timeToHitVerticalWall();
        double dtY = a.timeToHitHorizontalWall();
        if (t + dtX <= limit) 
        {
            pq.insert(new CollisionEvent(t + dtX, a, null));
        }
        if (t + dtY <= limit)
        {
            pq.insert(new CollisionEvent(t + dtY, null, a));
        } 
    }

    // redraw all particles
    private void redraw(double limit) {
        this.canvas.clear();
        for (int i = 0; i < particles.length; i++) {
            particles[i].draw(canvas);
        }
        this.canvas.show();

        if (t < limit) {
            pq.insert(new RedrawEvent(t + 1.0 / HZ));
        }
    }

    // infects particle S, adding its future recovery event and possible spread events to the queue
    private void infect(Particle S, double limit)
    {
        S.infected = true;
        S.inside_spread_zone = false;

        for (Particle particle : particles)
        {
            if (particle != S && particle.susceptible()) predict_zone_entry(particle, S, limit);
        }
        pq.insert(new RecoveryEvent(t + S.recovery_time, S));

        SIR_atual = new SIRCounter(SIR_atual.S - 1, SIR_atual.I + 1, SIR_atual.R);
    }

      
    /**
     * Simulates the system of particles for the specified amount of time.
     *
     * @param  limit the amount of time
     */
    public void simulate(double limit)
    {
        // initialize PQ with collision events and redraw event
        pq = new MinPQ<Event>();
        for (int i = 0; i < particles.length; i++) 
        {
            predict(particles[i], limit);
        }
        pq.insert(new RedrawEvent(0));        // redraw event

        SIR_atual = new SIRCounter(particles.length, 0, 0); // initialize SIR counter

        ST<Double, SIRCounter> tempos = new ST<Double, SIRCounter>(); // initialize graph data storage

        double scale = 1; 
        // graph is always from x = 0 to x = scale
        // scale grows as necessary to accommodate all points
        this.graph.setXscale(0, scale);
        this.graph.setYscale(0, particles.length);

        // infect patient zero
        infect(particles[0], limit);
        pq.insert(new UpdateGraphEvent(0));

        double last_print = 0;

        // the main event-driven simulation loop
        while (!pq.isEmpty()) { 
            if (t - last_print > 5)
            {
                last_print = t;
                if (SIR_atual.R > 0)
                {
                    double r = 0;
                    int i_seen = 0;
                    for (Particle particle : particles)
                    {
                        if (particle.recovered)
                        {
                            r += (double)particle.infected_count;
                            i_seen++;
                        }
                    }
                    //StdOut.println("r = " + r/i_seen);
                }
            }

            // get impending event, discard if invalidated
            Event event = pq.delMin();

            // categorize impending event
            // possible event type are Collision, EnterSpreadZone, InsideSpreadZone, Recovery, UpdateGraph and Redraw

            // handles collisions
            if (event instanceof CollisionEvent)
            {
                CollisionEvent e = (CollisionEvent)event;
                if (!e.isValid()) continue;

                Particle a = e.a;
                Particle b = e.b;

                // physical collision, so update positions, and then simulation clock
                for (int i = 0; i < particles.length; i++)
                    particles[i].move(e.time - t);
                t = e.time;

                // process event
                if      (a != null && b != null) a.bounceOff(b);              // particle-particle collision
                else if (a != null && b == null) a.bounceOffVerticalWall();   // particle-wall collision
                else if (a == null && b != null) b.bounceOffHorizontalWall(); // particle-wall collision

                // update the priority queue with new collisions involving a or b
                predict(a, limit);
                predict(b, limit);
            }

            // handles a healthy particle S entering I's spread zone 
            if (event instanceof EnterSpreadZoneEvent)
            {
                EnterSpreadZoneEvent e = (EnterSpreadZoneEvent)event;
                if (!e.isValid()) continue;

                Particle S = e.S;
                Particle I = e.I;

                S.inside_spread_zone = true;

                pq.insert(new InsideSpreadZoneEvent(t, S, I));
            }

            // triggers repeatedly while a healthy particle S is inside I's spread zone
            // possibly infects S and detects if it is still inside the spread zone 
            if (event instanceof InsideSpreadZoneEvent)
            {
                InsideSpreadZoneEvent e = (InsideSpreadZoneEvent)event;
                if (!e.isValid()) continue;

                Particle S = e.S;
                Particle I = e.I;

                // advances time to check if particle is still in spread zone
                for (int i = 0; i < particles.length; i++)
                    particles[i].move(e.time - t);
                t = e.time;

                double distance = Math.pow(S.rx - I.rx, 2) + Math.pow(S.ry - I.ry, 2);
                double max_distance = Math.pow(S.radius + I.spread_radius, 2);
                double min_distance = Math.pow(S.radius + I.radius, 2);

                // very small margin helps with small deviations
                if (max_distance - distance < -0.00009 || S.infected || I.recovered)
                {
                    S.inside_spread_zone = false;
                }
                else
                {
                    double i = -Math.log(INFECTIVITY)/(max_distance-min_distance);

                    double probability = BASE_INFECTIVITY * Math.exp(-i * (distance - min_distance)) / HZ_SPREAD_CHECK;
                    if (StdRandom.uniform(0.0, 1.0) < probability)
                    {
                        infect(S, limit);
                        I.infected_count++;
                        pq.insert(new UpdateGraphEvent(e.time));
                    }
                    else
                    {
                        pq.insert(new InsideSpreadZoneEvent(t + 1.0/HZ_SPREAD_CHECK, S, I));
                    }
                }
            }

            if (event instanceof RecoveryEvent)
            {
                RecoveryEvent e = (RecoveryEvent)event;

                e.I.recover();
                SIR_atual = new SIRCounter(SIR_atual.S, SIR_atual.I - 1, SIR_atual.R + 1);
                pq.insert(new UpdateGraphEvent(e.time));
            }
            
            if (event instanceof UpdateGraphEvent)
            {
                tempos.put(event.time, SIR_atual);
                if(t >= scale)
                {
                    while (t >= scale)
                    {
                        scale *= 1.1;
                    }
                    
                    this.graph.clear();
                    this.graph.setXscale(0, scale);
                    double last_time = 0;
                    for (double tempo : tempos.keys())
                    {
                        SIRCounter proximo = tempos.get(tempo);
                        SIRCounter ultimo = tempos.get(last_time);
                        this.graph.setPenColor(Draw.BLACK);
                        this.graph.line(last_time, ultimo.S, tempo, proximo.S);
                        this.graph.setPenColor(Draw.RED);
                        this.graph.line(last_time, ultimo.I, tempo, proximo.I);
                        this.graph.setPenColor(Draw.GRAY);
                        this.graph.line(last_time, ultimo.R, tempo, proximo.R);
                        last_time = tempo;
                    }
                }
                if (t > 0.000001)
                {
                    double ultimo_tempo = tempos.floor(t - 0.000001);
                    SIRCounter SIR_ultimo = tempos.get(ultimo_tempo);

                    this.graph.setPenColor(Draw.BLACK);
                    this.graph.line(ultimo_tempo, SIR_ultimo.S, t, SIR_atual.S);
                    this.graph.setPenColor(Draw.RED);
                    this.graph.line(ultimo_tempo, SIR_ultimo.I, t, SIR_atual.I);
                    this.graph.setPenColor(Draw.GRAY);
                    this.graph.line(ultimo_tempo, SIR_ultimo.R, t, SIR_atual.R);
                }
                this.graph.show();
            }


            if (event instanceof RedrawEvent)
            {
                RedrawEvent e = (RedrawEvent)event;
                for (int i = 0; i < particles.length; i++)
                    particles[i].move(e.time - t);
                t = e.time;
                redraw(limit);
            }
        }
    }

    abstract class Event implements Comparable<Event>
    {
        protected double time;

        public int compareTo(Event that) 
        {
            return Double.compare(this.time, that.time);
        }
    }

    class CollisionEvent extends Event
    {
        private final Particle a, b;
        private final long countA, countB;

        public CollisionEvent(double t, Particle a, Particle b) 
        {
            this.time = t;
            this.a    = a;
            this.b    = b;
            if (a != null) countA = a.count();
            else           countA = -1;
            if (b != null) countB = b.count();
            else           countB = -1;
        }

        public boolean isValid() {
            if (a != null && a.count() != countA) return false;
            if (b != null && b.count() != countB) return false;
            return true;
        }
    }

    class EnterSpreadZoneEvent extends Event
    {
        private final Particle S, I;
        private final long countA, countB;

        public EnterSpreadZoneEvent(double t, Particle S, Particle I) 
        {
            this.time = t;
            this.S    = S;
            this.I    = I;
            countA = S.count();
            countB = I.count();
        }

        public boolean isValid() {
            if (S.infected && I.susceptible()) return false;
            if (S != null && S.count() != countA) return false;
            if (I != null && I.count() != countB) return false;
            return true;
        }
    }

    class InsideSpreadZoneEvent extends Event
    {
        private final Particle S, I;

        public InsideSpreadZoneEvent (double t, Particle S, Particle I)
        {
            this.time = t;
            this.S = S;
            this.I = I;
        }

        public boolean isValid()
        {
            if (!S.inside_spread_zone)
            {
                return false;
            }
            return true;
        }
    }

    class RecoveryEvent extends Event
    {
        private final Particle I;

        public RecoveryEvent(double t, Particle I)
        {
            this.time = t;
            this.I = I;
        }
    }

    class UpdateGraphEvent extends Event
    {
        public UpdateGraphEvent(double t)
        {
            this.time = t;
        }
    }

    class RedrawEvent extends Event
    {
        public RedrawEvent(double t)
        {
            this.time = t;
        }
    }



    /**
     * Unit tests the {@code InfectionSim} data type.
     * Reads in the particle collision system from a standard input
     * (or generates {@code N} random particles if a command-line integer
     * is specified); simulates the system.
     
     * @param args the command-line arguments
     */
    public static void main(String[] args) {

        

        // the array of particles
        Particle[] particles;

        // create n random particles
        if (args.length == 1) {
            int n = Integer.parseInt(args[0]);
            particles = new Particle[n];
            for (int i = 0; i < n; i++)
                particles[i] = new Particle();
        }

        // or read from standard input
        else {
            int n = StdIn.readInt();
            particles = new Particle[n];
            StdRandom.setSeed(StdIn.readInt());

            INFECTIVITY = StdIn.readDouble();
            BASE_INFECTIVITY = StdIn.readDouble();

            for (int i = 0; i < n; i++) {
                double rx     = StdIn.readDouble();
                double ry     = StdIn.readDouble();
                double vx     = StdIn.readDouble();
                double vy     = StdIn.readDouble();
                double radius = StdIn.readDouble();
                double mass   = StdIn.readDouble();
                double recovery_time = StdIn.readDouble();
                double spread_radius = StdIn.readDouble();
                particles[i] = new Particle(rx, ry, vx, vy, radius, mass, Draw.BLACK, spread_radius, recovery_time);
            }
        }

        // create collision system and simulate
        InfectionSim system = new InfectionSim(particles);
        system.canvas.setCanvasSize(600, 600);

        // enable double buffering
        system.canvas.enableDoubleBuffering();

        system.simulate(10000);
    }
      
}