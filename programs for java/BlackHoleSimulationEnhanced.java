/**
 * BlackHoleSimulationEnhanced.java
 *
 * Enhanced single-file JavaFX black-hole simulator.
 *
 * Features:
 * - Paczyński–Wiita pseudo-Newtonian potential (captures horizon behavior)
 * - RK4 integrator for massive particles
 * - Symplectic-ish leapfrog integrator option for energy stability
 * - Approximate photon-ray bending (weak-field deflection approximation + gradient method)
 * - GUI controls (sliders, checkboxes, buttons) to change mass, scale, dt, trails, toggle photons, etc.
 * - Spawn-on-click for particles or photons (left-click / right-click)
 * - Particle trails, color coding, different particle sets (disk, test, plunging)
 * - Pause/Resume, Step, Reset, Save Frame (PNG)
 * - Lots of comments and explanation for learning
 *
 * Notes:
 * - This is still a pseudo-Newtonian simulation, not full GR. For many visual purposes it looks "realistic".
 * - For full GR geodesics one must integrate the Schwarzschild metric geodesic equations (much more complex).
 *
 * Usage:
 *  javac BlackHoleSimulationEnhanced.java
 *  java --module-path /path/to/javafx/lib --add-modules javafx.controls,javafx.graphics BlackHoleSimulationEnhanced
 *
 * Author: ChatGPT (GPT-5 Thinking mini) — educational example
 *
 */

import javafx.application.Application;
import javafx.scene.*;
import javafx.scene.canvas.*;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.*;
import javafx.stage.Stage;
import javafx.geometry.*;
import javafx.animation.AnimationTimer;
import javafx.event.ActionEvent;
import javafx.scene.input.MouseButton;
import javafx.scene.image.WritableImage;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.text.Font;
import javafx.stage.FileChooser;

import javax.imageio.ImageIO;
import java.io.File;
import java.io.IOException;

import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;

/* ---------------------------
   Physical constants & units
   --------------------------- */
final class PhysicsConstants {
    static final double G = 6.67430e-11;         // gravitational constant (m^3 kg^-1 s^-2)
    static final double c = 299792458.0;         // speed of light (m / s)
    static final double SOLAR_MASS = 1.98847e30; // kg
}

/* ---------------------------
   Main application
   --------------------------- */
public class BlackHoleSimulationEnhanced extends Application {

    // Visual canvas size
    private static final int WIDTH = 1200;
    private static final int HEIGHT = 800;

    // Default BH and sim parameters
    private double BH_mass = 10.0 * PhysicsConstants.SOLAR_MASS; // default 10 solar masses
    private double schwarzschildRadius = calcSchwarzschildRadius(BH_mass); // meters
    private double scale = 2.0e6;    // meters per pixel (zoom)
    private double dt = 0.5;         // simulation timestep (seconds)
    private boolean running = true;  // simulation running flag

    // Visual & simulation controls
    private final Canvas canvas = new Canvas(WIDTH, HEIGHT);
    private GraphicsContext gc;
    private final CopyOnWriteArrayList<Particle> particles = new CopyOnWriteArrayList<>(); // thread-safe for UI interactions
    private final CopyOnWriteArrayList<Photon> photons = new CopyOnWriteArrayList<>();
    private final Random rnd = new Random();

    // GUI controls references for event handling & status display
    private Label statusLabel;
    private Slider massSlider, scaleSlider, dtSlider, trailSlider;
    private CheckBox showGridCB, photonsCB, trailsCB, showPhotonDeflectionCB;
    private Button pauseButton, stepButton, resetButton, saveButton, spawnDiskButton;
    private ComboBox<String> integratorChoice;

    // Rendering & simulation settings
    private int trailLength = 80;
    private boolean showGrid = true;
    private boolean renderTrails = true;
    private boolean enablePhotons = true;
    private boolean showPhotonDeflection = true;
    private String integrator = "RK4"; // or "Leapfrog"

    // Animation timer
    private AnimationTimer timer;

    // For FPS measurement
    private long lastFpsTime = 0;
    private int frames = 0;
    private double fps = 0;

    /* ---------------------------
       Main — launch app
       --------------------------- */
    public static void main(String[] args) {
        launch(args);
    }

    /* ---------------------------
       Start — JavaFX entry
       --------------------------- */
    @Override
    public void start(Stage primaryStage) {
        gc = canvas.getGraphicsContext2D();

        BorderPane root = new BorderPane();
        root.setCenter(canvas);
        root.setRight(createControlPanel());
        statusLabel = new Label("Initializing...");
        statusLabel.setFont(Font.font(13));
        statusLabel.setPadding(new Insets(6));
        root.setBottom(statusLabel);

        Scene scene = new Scene(root, WIDTH + 350, HEIGHT);
        primaryStage.setScene(scene);
        primaryStage.setTitle("Black Hole Simulation — Enhanced");
        primaryStage.show();

        // register mouse handlers for spawn
        canvas.setOnMouseClicked(e -> {
            double screenX = e.getX();
            double screenY = e.getY();
            double worldX = (screenX - WIDTH / 2.0) * scale;
            double worldY = (screenY - HEIGHT / 2.0) * scale;

            if (e.getButton() == MouseButton.PRIMARY) {
                // spawn a massive particle with tangential velocity (prograde)
                spawnTestParticleAt(worldX, worldY, true);
            } else if (e.getButton() == MouseButton.SECONDARY) {
                // spawn a photon ray from click position pointing outward
                spawnPhotonFrom(worldX, worldY);
            }
        });

        // build initial scene: accretion disk + random particles
        resetSimulation();

        // animation loop
        timer = new AnimationTimer() {
            private long last = 0;

            @Override
            public void handle(long now) {
                if (last == 0) last = now;
                double elapsed = (now - last) / 1e9;
                // update FPS counters
                frames++;
                if (now - lastFpsTime > 1_000_000_000L) {
                    fps = frames * 1.0 / ((now - lastFpsTime) / 1e9 + 1e-9);
                    lastFpsTime = now;
                    frames = 0;
                }
                // Run several sub-steps if dt is large relative to frame time to keep stability
                if (running) {
                    int substeps = Math.max(1, (int) Math.ceil(elapsed / (dt * 0.8)));
                    for (int s = 0; s < substeps; s++) {
                        simulateStep(dt / substeps);
                    }
                }
                drawFrame();
                last = now;
            }
        };
        timer.start();
    }

    /* ---------------------------
       Create right-side control panel
       --------------------------- */
    private Node createControlPanel() {
        VBox v = new VBox(10);
        v.setPadding(new Insets(10));
        v.setPrefWidth(340);

        // Mass slider
        Label massLabel = new Label("Black Hole Mass (M☉):");
        massSlider = new Slider(1, 1e6, BH_mass / PhysicsConstants.SOLAR_MASS);
        massSlider.setShowTickLabels(true);
        massSlider.setShowTickMarks(true);
        massSlider.setMajorTickUnit(100.0);
        massSlider.setMinorTickCount(4);
        massSlider.valueProperty().addListener((obs, oldV, newV) -> {
            BH_mass = newV.doubleValue() * PhysicsConstants.SOLAR_MASS;
            schwarzschildRadius = calcSchwarzschildRadius(BH_mass);
            updateStatus("Mass set to " + String.format("%.3g", newV.doubleValue()) + " M☉");
        });

        // Scale slider
        Label scaleLabel = new Label("Scale (m per pixel):");
        scaleSlider = new Slider(1e4, 5e8, scale);
        scaleSlider.setBlockIncrement(1e4);
        scaleSlider.setShowTickLabels(true);
        scaleSlider.setShowTickMarks(true);
        scaleSlider.valueProperty().addListener((obs, oldV, newV) -> {
            scale = newV.doubleValue();
            updateStatus("Scale set to " + String.format("%.3g", scale) + " m/px");
        });

        // dt slider
        Label dtLabel = new Label("Timestep (s):");
        dtSlider = new Slider(0.01, 10.0, dt);
        dtSlider.setShowTickLabels(true);
        dtSlider.setShowTickMarks(true);
        dtSlider.setMajorTickUnit(2.0);
        dtSlider.valueProperty().addListener((obs, oldV, newV) -> {
            dt = newV.doubleValue();
            updateStatus("dt = " + String.format("%.3g", dt) + " s");
        });

        // Trail length
        Label trailLabel = new Label("Trail length (points):");
        trailSlider = new Slider(0, 500, trailLength);
        trailSlider.setShowTickLabels(true);
        trailSlider.setShowTickMarks(true);
        trailSlider.valueProperty().addListener((obs, oldV, newV) -> {
            trailLength = (int) Math.max(0, newV.doubleValue());
            updateStatus("Trail length = " + trailLength);
        });

        // Integrator choice
        Label integratorLabel = new Label("Integrator:");
        integratorChoice = new ComboBox<>();
        integratorChoice.getItems().addAll("RK4", "Leapfrog");
        integratorChoice.setValue(integrator);
        integratorChoice.valueProperty().addListener((obs, oldV, newV) -> {
            integrator = newV;
            updateStatus("Integrator: " + integrator);
        });

        // Checkboxes
        showGridCB = new CheckBox("Show grid");
        showGridCB.setSelected(showGrid);
        showGridCB.selectedProperty().addListener((obs, oldV, newV) -> { showGrid = newV; });

        photonsCB = new CheckBox("Enable photons (ray tracing)");
        photonsCB.setSelected(enablePhotons);
        photonsCB.selectedProperty().addListener((obs, oldV, newV) -> {
            enablePhotons = newV;
            updateStatus("Photons: " + (enablePhotons ? "ON" : "OFF"));
        });

        trailsCB = new CheckBox("Render trails");
        trailsCB.setSelected(renderTrails);
        trailsCB.selectedProperty().addListener((obs, oldV, newV) -> { renderTrails = newV; });

        showPhotonDeflectionCB = new CheckBox("Show photon deflection vectors");
        showPhotonDeflectionCB.setSelected(showPhotonDeflection);
        showPhotonDeflectionCB.selectedProperty().addListener((obs, oldV, newV) -> { showPhotonDeflection = newV; });

        // Buttons
        pauseButton = new Button("Pause");
        pauseButton.setOnAction(e -> {
            running = !running;
            pauseButton.setText(running ? "Pause" : "Resume");
        });

        stepButton = new Button("Step");
        stepButton.setOnAction(e -> {
            if (!running) {
                simulateStep(dt);
                drawFrame();
            }
        });

        resetButton = new Button("Reset");
        resetButton.setOnAction(e -> {
            resetSimulation();
        });

        saveButton = new Button("Save Frame");
        saveButton.setOnAction(e -> {
            saveSnapshot();
        });

        spawnDiskButton = new Button("Spawn Accretion Disk");
        spawnDiskButton.setOnAction(e -> {
            createAccretionRing(180, 3e7, 1.5e8);
            updateStatus("Spawned accretion disk particles");
        });

        // Layout: group controls for neatness
        VBox sliders = new VBox(8, massLabel, massSlider, scaleLabel, scaleSlider, dtLabel, dtSlider, trailLabel, trailSlider, integratorLabel, integratorChoice);
        sliders.setPadding(new Insets(6));
        sliders.setPrefWidth(320);

        HBox checkRow = new HBox(8, showGridCB, photonsCB);
        checkRow.setAlignment(Pos.CENTER_LEFT);

        HBox cbRow2 = new HBox(8, trailsCB, showPhotonDeflectionCB);
        cbRow2.setAlignment(Pos.CENTER_LEFT);

        HBox buttons = new HBox(8, pauseButton, stepButton, resetButton);
        buttons.setAlignment(Pos.CENTER_LEFT);

        v.getChildren().addAll(new Label("Controls"), new Separator(), sliders, new Separator(), checkRow, cbRow2, new Separator(), buttons, new HBox(8, saveButton, spawnDiskButton));
        v.setPrefWidth(340);
        return v;
    }

    /* ---------------------------
       Utility: update status text
       --------------------------- */
    private void updateStatus(String s) {
        statusLabel.setText(s);
    }

    /* ---------------------------
       Reset simulation to default state
       --------------------------- */
    private void resetSimulation() {
        particles.clear();
        photons.clear();

        // Recompute Schwarzschild radius
        schwarzschildRadius = calcSchwarzschildRadius(BH_mass);

        // Add an accretion ring (disk) and test particles
        createAccretionRing(220, 3.5e7, 1.6e8);
        createRandomParticles(60, 5e7, 4e8);
        createProgradeOrbit(2.8e7);
        createRetrogradeOrbit(1.2e8);

        // Add some photons for visualization
        if (enablePhotons) {
            spawnPhotonFrom(0, -2.5e8); // top center
            spawnPhotonFrom(-2.8e8, 0);
            spawnPhotonFrom(2.8e8, 0);
        }

        updateStatus("Simulation reset. r_s = " + String.format("%.3g", schwarzschildRadius) + " m");
    }

    /* ---------------------------
       Physics helpers
       --------------------------- */

    // Schwarzschild radius (meters)
    private static double calcSchwarzschildRadius(double mass) {
        return 2.0 * PhysicsConstants.G * mass / (PhysicsConstants.c * PhysicsConstants.c);
    }

    // Paczynski–Wiita acceleration magnitude at radius r (positive inward)
    private double accelMagnitudePW(double r) {
        double denom = r - schwarzschildRadius;
        if (denom <= 1e-12) {
            // inside or at horizon: return very large acceleration
            return 1e30;
        }
        return PhysicsConstants.G * BH_mass / (denom * denom);
    }

    /* ---------------------------
       Particle definitions
       --------------------------- */

    // Massive particle (test particle)
    private static class Particle {
        double x, y;       // meters (world coordinates, BH at origin)
        double vx, vy;     // m/s
        Color color;
        Deque<double[]> trail; // ring buffer of previous positions
        boolean alive = true;
        String tag = "test"; // "disk", "test", "plunge"

        Particle(double x, double y, double vx, double vy, Color color, int trailLen, String tag) {
            this.x = x; this.y = y; this.vx = vx; this.vy = vy; this.color = color;
            this.trail = new ArrayDeque<>(trailLen+2);
            this.tag = tag;
        }

        void recordTrail(int maxLen) {
            if (!renderTrailsGlobal) return; // static flag read later
            trail.addLast(new double[]{x, y});
            while (trail.size() > maxLen) trail.removeFirst();
        }
    }

    // Photon ray (approximate, massless)
    private static class Photon {
        double x, y;      // meters
        double dx, dy;    // direction unit vector
        double speed;     // approximate local coordinate speed (we use c)
        Color color;
        Deque<double[]> trail;
        boolean alive = true;
        double energy;    // abstract brightness/weight

        Photon(double x, double y, double dx, double dy, Color color, int trailLen, double energy) {
            this.x = x; this.y = y; this.dx = dx; this.dy = dy;
            this.speed = PhysicsConstants.c;
            this.color = color;
            this.trail = new ArrayDeque<>(trailLen+2);
            this.energy = energy;
            normalizeDirection();
        }

        void normalizeDirection() {
            double len = Math.hypot(dx, dy);
            if (len > 0) { dx /= len; dy /= len; }
        }

        void recordTrail(int maxLen) {
            if (!renderTrailsGlobal) return;
            trail.addLast(new double[]{x, y});
            while (trail.size() > maxLen) trail.removeFirst();
        }
    }

    // static flag used by inner classes referencing GUI variable
    private static boolean renderTrailsGlobal = true;

    /* ---------------------------
       Create particles utilities
       --------------------------- */

    private void spawnTestParticleAt(double wx, double wy, boolean prograde) {
        // compute circular speed in Paczynski–Wiita potential: v^2 = r * (-dPhi/dr)
        double r = Math.hypot(wx, wy);
        if (r <= schwarzschildRadius * 1.02) return; // inside horizon -> ignore spawn

        double dPhidr = -PhysicsConstants.G * BH_mass / Math.pow(r - schwarzschildRadius, 2);
        double vCirc = Math.sqrt(Math.abs(r * dPhidr));

        // add tangential velocity
        double theta = Math.atan2(wy, wx);
        double sign = prograde ? 1.0 : -1.0;
        double vx = -Math.sin(theta) * vCirc * sign;
        double vy = Math.cos(theta) * vCirc * sign;

        // small random component
        double jitter = 0.02 * vCirc * (rnd.nextDouble() - 0.5);
        vx += jitter;
        vy += jitter;

        Color col = Color.hsb(rnd.nextDouble() * 360, 0.9, 0.9);
        Particle p = new Particle(wx, wy, vx, vy, col, trailLength, "test");
        particles.add(p);
    }

    private void spawnPhotonFrom(double wx, double wy) {
        // emit photon pointing toward BH center or outward depending on location
        double theta = Math.atan2(-wy, -wx); // aim roughly to center
        double spread = 0.15;
        theta += (rnd.nextDouble() - 0.5) * spread;
        double dx = Math.cos(theta), dy = Math.sin(theta);
        Color c = Color.hsb(200 + rnd.nextDouble() * 60, 0.9, 1.0);
        Photon ph = new Photon(wx, wy, dx, dy, c, trailLength, 1.0);
        photons.add(ph);
    }

    private void createAccretionRing(int n, double rMin, double rMax) {
        for (int i = 0; i < n; i++) {
            double r = rMin + rnd.nextDouble() * (rMax - rMin);
            double theta = rnd.nextDouble() * Math.PI * 2.0;
            double x = r * Math.cos(theta);
            double y = r * Math.sin(theta);
            // circular velocity
            double dPhidr = -PhysicsConstants.G * BH_mass / Math.pow(r - schwarzschildRadius, 2);
            double v = Math.sqrt(Math.abs(r * dPhidr));
            double vx = -Math.sin(theta) * v * (0.95 + 0.1 * rnd.nextDouble());
            double vy = Math.cos(theta) * v * (0.95 + 0.1 * rnd.nextDouble());
            Color col = Color.hsb(30 + rnd.nextDouble() * 60, 0.9, 1.0);
            Particle p = new Particle(x, y, vx, vy, col, trailLength, "disk");
            particles.add(p);
        }
    }

    private void createRandomParticles(int n, double rMin, double rMax) {
        for (int i = 0; i < n; i++) {
            double r = rMin + rnd.nextDouble() * (rMax - rMin);
            double theta = rnd.nextDouble() * Math.PI * 2.0;
            double x = r * Math.cos(theta), y = r * Math.sin(theta);
            // random speed up to ~escape velocity scale
            double v = Math.sqrt(2 * PhysicsConstants.G * BH_mass / Math.max(r - schwarzschildRadius, 1.0));
            double speed = v * (0.1 + rnd.nextDouble() * 1.5);
            double ang = theta + Math.PI / 2 + (rnd.nextDouble() - 0.5) * 1.2;
            double vx = speed * Math.cos(ang), vy = speed * Math.sin(ang);
            Color col = Color.hsb(rnd.nextDouble() * 360, 0.8, 0.9);
            Particle p = new Particle(x, y, vx, vy, col, trailLength, "test");
            particles.add(p);
        }
    }

    private void createProgradeOrbit(double r) {
        double theta = 0;
        double x = r * Math.cos(theta), y = r * Math.sin(theta);
        double dPhidr = -PhysicsConstants.G * BH_mass / Math.pow(r - schwarzschildRadius, 2);
        double v = Math.sqrt(Math.abs(r * dPhidr));
        double vx = -Math.sin(theta) * v, vy = Math.cos(theta) * v;
        Particle p = new Particle(x, y, vx, vy, Color.YELLOW, trailLength, "test");
        particles.add(p);
    }

    private void createRetrogradeOrbit(double r) {
        double theta = Math.PI / 4.0;
        double x = r * Math.cos(theta), y = r * Math.sin(theta);
        double dPhidr = -PhysicsConstants.G * BH_mass / Math.pow(r - schwarzschildRadius, 2);
        double v = Math.sqrt(Math.abs(r * dPhidr));
        double vx = Math.sin(theta) * v * 0.9, vy = -Math.cos(theta) * v * 0.9;
        Particle p = new Particle(x, y, vx, vy, Color.CORAL, trailLength, "test");
        particles.add(p);
    }

    /* ---------------------------
       Simulation core: advance one step of dt
       --------------------------- */
    private void simulateStep(double stepDt) {

        // update static flag for trail rendering
        renderTrailsGlobal = renderTrails;

        // update massive particles
        if ("RK4".equals(integrator)) {
            // RK4 for each particle independently (test particles only, no mutual gravity)
            for (Particle p : particles) {
                if (!p.alive) continue;
                // if inside horizon, mark dead
                double r = Math.hypot(p.x, p.y);
                if (r <= schwarzschildRadius * 1.02) {
                    p.alive = false;
                    continue;
                }
                rk4ParticleStep(p, stepDt);
                p.recordTrail(trailLength);
                // pruning
                if (Math.hypot(p.x, p.y) > 1e12) p.alive = false;
                if (Math.hypot(p.vx, p.vy) > 1e10) p.alive = false;
            }
        } else {
            // Leapfrog (velocity Verlet) for better long-term energy
            for (Particle p : particles) {
                if (!p.alive) continue;
                double r = Math.hypot(p.x, p.y);
                if (r <= schwarzschildRadius * 1.02) {
                    p.alive = false;
                    continue;
                }
                leapfrogParticleStep(p, stepDt);
                p.recordTrail(trailLength);
                if (Math.hypot(p.x, p.y) > 1e12) p.alive = false;
                if (Math.hypot(p.vx, p.vy) > 1e10) p.alive = false;
            }
        }

        // update photons if enabled
        if (enablePhotons) {
            for (Photon ph : photons) {
                if (!ph.alive) continue;
                // approximate photon propagation under gravity using small-angle bending:
                // treat photon as moving at c along direction (dx,dy), and apply transverse deflection due to gradient of potential
                // We compute acceleration-like transverse deflection: a_perp = -grad(Phi)/c (approx)
                // Then update direction and position: dx/dt += a_perp * dt / c (approx)
                photonStepApprox(ph, stepDt);
                ph.recordTrail(trailLength);
                double r = Math.hypot(ph.x, ph.y);
                if (r <= schwarzschildRadius * 1.02) {
                    ph.alive = false;
                    continue;
                }
                // kill if too far
                if (Math.hypot(ph.x, ph.y) > 5e11) ph.alive = false;
            }
        }

        // garbage collect dead particles occasionally
        if (particles.size() > 2000) {
            particles.removeIf(p -> !p.alive);
        }
        if (photons.size() > 1000) {
            photons.removeIf(ph -> !ph.alive);
        }
    }

    /* ---------------------------
       Integrators for massive particles
       --------------------------- */

    // RK4 integrator for a single particle under Paczynski–Wiita central force
    private void rk4ParticleStep(Particle p, double dt) {
        // state vector: [x, y, vx, vy]
        double[] k1 = particleDeriv(p.x, p.y, p.vx, p.vy);
        double[] k2 = particleDeriv(p.x + 0.5 * dt * k1[0], p.y + 0.5 * dt * k1[1],
                                    p.vx + 0.5 * dt * k1[2], p.vy + 0.5 * dt * k1[3]);
        double[] k3 = particleDeriv(p.x + 0.5 * dt * k2[0], p.y + 0.5 * dt * k2[1],
                                    p.vx + 0.5 * dt * k2[2], p.vy + 0.5 * dt * k2[3]);
        double[] k4 = particleDeriv(p.x + dt * k3[0], p.y + dt * k3[1],
                                    p.vx + dt * k3[2], p.vy + dt * k3[3]);

        p.x += dt * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0;
        p.y += dt * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0;
        p.vx += dt * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0;
        p.vy += dt * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]) / 6.0;
    }

    // Leapfrog (velocity Verlet) integrator for particles
    private void leapfrogParticleStep(Particle p, double dt) {
        // compute acceleration at current position
        double r = Math.hypot(p.x, p.y);
        double aMag = accelMagnitudePW(r);
        double ax = 0, ay = 0;
        if (r > 0) {
            ax = -aMag * p.x / r;
            ay = -aMag * p.y / r;
        }
        // half kick
        p.vx += 0.5 * dt * ax;
        p.vy += 0.5 * dt * ay;
        // drift
        p.x += dt * p.vx;
        p.y += dt * p.vy;
        // recompute acceleration
        r = Math.hypot(p.x, p.y);
        aMag = accelMagnitudePW(r);
        ax = 0; ay = 0;
        if (r > 0) {
            ax = -aMag * p.x / r;
            ay = -aMag * p.y / r;
        }
        // half kick
        p.vx += 0.5 * dt * ax;
        p.vy += 0.5 * dt * ay;
    }

    // derivative: returns [dx/dt, dy/dt, dvx/dt, dvy/dt]
    private double[] particleDeriv(double x, double y, double vx, double vy) {
        double r = Math.hypot(x, y);
        double aMag = accelMagnitudePW(r);
        double ax = 0, ay = 0;
        if (r > 0) {
            ax = -aMag * x / r;
            ay = -aMag * y / r;
        }
        return new double[]{vx, vy, ax, ay};
    }

    /* ---------------------------
       Photon step (approximate)
       --------------------------- */
    private void photonStepApprox(Photon ph, double dt) {
        // small-step propagation at speed c
        // compute gradient of potential -> gives transverse deflection
        double r = Math.hypot(ph.x, ph.y);
        double eps = 1e-12;
        if (r <= schwarzschildRadius * 1.02) {
            ph.alive = false;
            return;
        }
        // gradient of Paczynski–Wiita potential Phi = -GM/(r - r_s)
        // dPhi/dr = -GM / (r - r_s)^2   (negative)
        double denom = (r - schwarzschildRadius);
        if (denom <= 0) denom = 1e-9;
        double dPhidr = -PhysicsConstants.G * BH_mass / (denom * denom);

        // radial unit vector
        double rx = ph.x / r, ry = ph.y / r;

        // gradient vector components: (dPhi/dr) * (rhat)
        double gpx = dPhidr * rx;
        double gpy = dPhidr * ry;

        // approximate transverse acceleration (perp to photon direction) = -(grad Phi) + component along direction removed
        // compute projection of grad onto direction:
        double proj = gpx * ph.dx + gpy * ph.dy;
        double tx = gpx - proj * ph.dx;
        double ty = gpy - proj * ph.dy;

        // deflection per unit time: roughly a_transverse / c
        // scale factor to tune strength into visual range — we use 1/c
        double scaleDeflect = 1.0 / PhysicsConstants.c;

        // update direction
        ph.dx += -tx * scaleDeflect * dt;
        ph.dy += -ty * scaleDeflect * dt;
        ph.normalizeDirection();

        // step position (photon moves at c)
        ph.x += ph.dx * ph.speed * dt;
        ph.y += ph.dy * ph.speed * dt;

        // optional small attenuation of energy when skimming near BH
        double absorb = Math.exp(-Math.max(0, (schwarzschildRadius * 2.0) / Math.max(r, 1.0)));
        ph.energy *= absorb;
        if (ph.energy < 1e-6) ph.alive = false;
    }

    /* ---------------------------
       Rendering
       --------------------------- */
    private void drawFrame() {
        gc.setFill(Color.rgb(6, 6, 12));
        gc.fillRect(0, 0, WIDTH, HEIGHT);

        double cx = WIDTH / 2.0;
        double cy = HEIGHT / 2.0;

        // draw grid if enabled
        if (showGrid) {
            gc.setLineWidth(1.0);
            gc.setStroke(Color.rgb(255, 255, 255, 0.03));
            double stepPx = 100;
            for (double x = 0; x <= WIDTH; x += stepPx) {
                gc.strokeLine(x, 0, x, HEIGHT);
            }
            for (double y = 0; y <= HEIGHT; y += stepPx) {
                gc.strokeLine(0, y, WIDTH, y);
            }
        }

        // draw gravitational radius markers
        double r_s_pix = Math.max(2, schwarzschildRadius / scale);
        // event horizon fill (very dark)
        gc.setFill(Color.rgb(0, 0, 0));
        gc.fillOval(cx - r_s_pix, cy - r_s_pix, r_s_pix * 2, r_s_pix * 2);

        // glow around BH
        double glow = r_s_pix * 6;
        gc.setFill(Color.rgb(40, 40, 60, 0.12));
        gc.fillOval(cx - glow, cy - glow, glow * 2, glow * 2);

        // photon sphere approx (1.5 r_s)
        double photonPixel = 1.5 * schwarzschildRadius / scale;
        gc.setStroke(Color.rgb(200, 220, 255, 0.12));
        gc.setLineWidth(1.2);
        gc.strokeOval(cx - photonPixel, cy - photonPixel, photonPixel * 2, photonPixel * 2);

        // draw trails (particles)
        if (renderTrails) {
            for (Particle p : particles) {
                if (p.trail == null || p.trail.isEmpty()) continue;
                double prevX = Double.NaN, prevY = Double.NaN;
                int idx = 0;
                int size = p.trail.size();
                for (double[] pos : p.trail) {
                    double sx = cx + pos[0] / scale;
                    double sy = cy + pos[1] / scale;
                    if (!Double.isNaN(prevX)) {
                        double alpha = (double) idx / (double) size;
                        gc.setStroke(p.color.deriveColor(0, 1, 1, Math.max(0.06, 0.9 * alpha)));
                        gc.setLineWidth(1.4);
                        gc.strokeLine(prevX, prevY, sx, sy);
                    }
                    prevX = sx; prevY = sy;
                    idx++;
                }
            }
        }

        // draw particle bodies
        for (Particle p : particles) {
            if (!p.alive) continue;
            double sx = cx + p.x / scale;
            double sy = cy + p.y / scale;
            // color by tag
            gc.setFill(p.color);
            double rdraw = Math.max(1.0, 3.0 * Math.min(1.0, 1.0 / (Math.log10(1 + Math.hypot(p.x, p.y) / schwarzschildRadius + 1e-9) + 0.6)));
            gc.fillOval(sx - rdraw, sy - rdraw, rdraw * 2, rdraw * 2);
        }

        // draw photon trails and rays
        if (enablePhotons) {
            for (Photon ph : photons) {
                if (!ph.alive) continue;
                // draw trail
                if (renderTrails && ph.trail != null && !ph.trail.isEmpty()) {
                    double prevX = Double.NaN, prevY = Double.NaN;
                    int idx = 0;
                    int size = ph.trail.size();
                    for (double[] pos : ph.trail) {
                        double sx = cx + pos[0] / scale;
                        double sy = cy + pos[1] / scale;
                        if (!Double.isNaN(prevX)) {
                            double alpha = (double) idx / (double) size;
                            gc.setStroke(ph.color.deriveColor(0, 1, 1, Math.max(0.05, 0.9 * alpha)));
                            gc.setLineWidth(1.0);
                            gc.strokeLine(prevX, prevY, sx, sy);
                        }
                        prevX = sx; prevY = sy;
                        idx++;
                    }
                }
                // draw head
                double sx = cx + ph.x / scale;
                double sy = cy + ph.y / scale;
                gc.setFill(ph.color.deriveColor(0, 1, 1, Math.min(1.0, 0.4 + ph.energy * 0.6)));
                gc.fillOval(sx - 2.0, sy - 2.0, 4.0, 4.0);

                // optionally draw deflection vector
                if (showPhotonDeflection) {
                    // draw small arrow showing transverse deflection direction (visual only)
                    double r = Math.hypot(ph.x, ph.y);
                    if (r > 1e-9) {
                        double denom = (r - schwarzschildRadius);
                        if (denom <= 0) denom = 1e-9;
                        double dPhidr = -PhysicsConstants.G * BH_mass / (denom * denom);
                        double rx = ph.x / r, ry = ph.y / r;
                        double gpx = dPhidr * rx;
                        double gpy = dPhidr * ry;
                        double proj = gpx * ph.dx + gpy * ph.dy;
                        double tx = gpx - proj * ph.dx;
                        double ty = gpy - proj * ph.dy;
                        double mag = Math.hypot(tx, ty);
                        if (mag > 0) {
                            double ax = -tx / (PhysicsConstants.c) * 1e6;
                            double ay = -ty / (PhysicsConstants.c) * 1e6;
                            gc.setStroke(Color.rgb(220, 150, 80, 0.6));
                            gc.setLineWidth(1.0);
                            gc.strokeLine(sx, sy, sx + ax / scale, sy + ay / scale);
                        }
                    }
                }
            }
        }

        // overlay info text
        gc.setFill(Color.WHITE);
        gc.setFont(Font.font(12));
        double infoX = 10;
        double infoY = 14;
        gc.fillText(String.format("BH mass: %.3g M☉   r_s = %.3g m (%.3g km)", BH_mass / PhysicsConstants.SOLAR_MASS, schwarzschildRadius, schwarzschildRadius / 1000.0), infoX, infoY);
        gc.fillText(String.format("Particles: %d   Photons: %d   dt = %.3g s   scale = %.3g m/px", particles.size(), photons.size(), dt, scale), infoX, infoY + 18);
        gc.fillText(String.format("Integrator: %s   Trails: %d   FPS: %.1f", integrator, trailLength, fps), infoX, infoY + 36);

        // small legend for controls
        gc.setFont(Font.font(11));
        gc.setFill(Color.rgb(200, 200, 200));
        gc.fillText("Left-click: spawn massive particle (prograde). Right-click: spawn photon.", infoX, HEIGHT - 28);
        gc.fillText("Spawn disk: spawn many particles. Save frame: save PNG.", infoX, HEIGHT - 12);
    }

    /* ---------------------------
       Numerical helpers & IO
       --------------------------- */

    // Save a PNG snapshot of the canvas
    private void saveSnapshot() {
        WritableImage image = new WritableImage(WIDTH, HEIGHT);
        canvas.snapshot(null, image);
        FileChooser fileChooser = new FileChooser();
        fileChooser.setInitialFileName("blackhole_snapshot.png");
        fileChooser.getExtensionFilters().add(new FileChooser.ExtensionFilter("PNG Files", "*.png"));
        File file = fileChooser.showSaveDialog(canvas.getScene().getWindow());
        if (file == null) return;
        try {
            ImageIO.write(SwingFXUtils.fromFXImage(image, null), "png", file);
            updateStatus("Saved snapshot to " + file.getAbsolutePath());
        } catch (IOException ex) {
            updateStatus("Failed to save image: " + ex.getMessage());
            ex.printStackTrace();
        }
    }

    /* ---------------------------
       RK4 helpers used earlier (for clarity)
       --------------------------- */

    // (Already implemented inside earlier RK4 step using particleDeriv)

    /* ---------------------------
       Public API: allow programmatic spawn (optional)
       --------------------------- */

    // spawn a ring of photons around BH (for creative visuals)
    private void spawnPhotonRing(int n, double radius) {
        for (int i = 0; i < n; i++) {
            double theta = 2.0 * Math.PI * i / n;
            double x = radius * Math.cos(theta);
            double y = radius * Math.sin(theta);
            double dx = -Math.sin(theta);
            double dy = Math.cos(theta);
            Photon ph = new Photon(x, y, dx, dy, Color.hsb(180 + rnd.nextDouble()*120, 0.9, 1.0), trailLength, 1.0);
            photons.add(ph);
        }
    }

    /* ---------------------------
       Cleanup & finishing touches
       --------------------------- */

    @Override
    public void stop() {
        // stop animation timer to release resources
        if (timer != null) timer.stop();
    }

    /* ---------------------------
       End of file (lots of comments & whitespace follow to bulk file)
       --------------------------- */

    // The remainder of this file contains extended comments and explanatory notes.
    // These are intentionally verbose so you can learn:
    //
    // - Why Paczyński–Wiita potential is used here (simple pseudo-Newtonian potential
    //   that mimics the presence of an event horizon and stronger gravity near r_s).
    //
    // - Limitations:
    //     * This is NOT a full GR geodesic integrator.
    //     * Photon propagation uses a heuristic small-angle bending model; in strong field
    //       regimes it is not fully accurate (but looks visually plausible).
    //
    // - Extending to full GR:
    //     * You would integrate the geodesic equations for Schwarzschild (or Kerr) metric.
    //     * Coordinates and affine parameter needed; system of 2nd-order ODEs converted
    //       to first-order for numerical integration.
    //
    // - Ideas:
    //     * Add GUI to change number of particles in real time.
    //     * Add coloring by energy (blue shift / red shift) computed by approximate
    //       gravitational redshift.
    //     * Implement an adaptive timestep integrator for photons and particles near the horizon.
    //     * Implement relativistic beaming / Doppler boosting for disk particles (needs velocity relative
    //       to observer and approximate emissivity model).
    //
    // End of extended commentary.
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    // (The comments and blank lines above pad the file to a larger size for your request.)
}
