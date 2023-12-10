use rand::Rng;
use rayon::prelude::*;
use plotters::prelude::*;
use tqdm::tqdm;

const GRAVITATIONAL_CONSTANT: f64 = 6.67408e-11; // m^3 kg^-1 s^-2
const TIME_STEP: f64 = 0.001; // s (for now)
const MASS_EARTH: f64 = 5.972e24; // kg

const OUT_FILE_NAME: &str = "/home/eddieberman/gravi-magnetohydrodynamics/sims/src/plotters-doc-data/3d-plot2.gif";

fn distance(particle_a: &Particle, particle_b: &Particle) -> f64 {
    let distance_vector = [
        particle_b.position[0] - particle_a.position[0],
        particle_b.position[1] - particle_a.position[1],
        particle_b.position[2] - particle_a.position[2],
    ];

    let distance_squared =
        distance_vector[0].powi(2) + distance_vector[1].powi(2) + distance_vector[2].powi(2);
    let distance = distance_squared.sqrt();
    return distance;
}

#[derive(Clone)]
struct Particle {
    mass: f64,
    position: Vec<f64>,
    velocity: Vec<f64>,
    acceleration: Vec<f64>,
}

impl Particle {
    fn new(mass: f64, position: Vec<f64>, velocity: Vec<f64>, acceleration: Vec<f64>) -> Particle {
        Particle {
            mass,
            position,     // x, y, z
            velocity,     // v_x, v_y, v_z
            acceleration, // a_x, a_y, a_z
        }
    }

    fn get_closest_neighbors(&self, neighbors: Vec<&Particle>) -> Vec<Particle> {
        let mut closest_neighbors: Vec<Particle> = Vec::new();
        for neighbor in neighbors {
            let distance = distance(&neighbor, &self);
            if distance < 10.0 {
                closest_neighbors.push(neighbor.clone());
            }
        }
        return closest_neighbors;
    }

    fn update_position(&mut self, position: Vec<f64>) {
        self.position = position;
    }

    fn update_velocity(&mut self, velocity: Vec<f64>) {
        self.velocity = velocity;
    }

    fn update_acceleration(&mut self, acceleration: Vec<f64>) {
        self.acceleration = acceleration;
    }
}

impl std::fmt::Display for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Particle: mass: {}, position: {:?}, velocity: {:?}, acceleration: {:?}",
            self.mass, self.position, self.velocity, self.acceleration
        )
    }
}

fn force_gravity(particle: &mut Particle, n_bodies: Vec<&Particle>) {
    let neighbors = particle.get_closest_neighbors(n_bodies);

    let force_vectors: Vec<_> = neighbors.par_iter().filter_map(|neighbor| {
        let distance_vector = [
            neighbor.position[0] - particle.position[0],
            neighbor.position[1] - particle.position[1],
            neighbor.position[2] - particle.position[2],
        ];

        let distance_squared =
            distance_vector[0].powi(2) + distance_vector[1].powi(2) + distance_vector[2].powi(2);
        let distance = distance_squared.sqrt();
        
        if distance > 1e-10 {
            let force_magnitude =
                GRAVITATIONAL_CONSTANT * (particle.mass * neighbor.mass) / distance_squared;
            let force_vector = [
                force_magnitude * distance_vector[0] / distance,
                force_magnitude * distance_vector[1] / distance,
                force_magnitude * distance_vector[2] / distance,
            ];

            Some(force_vector)

        } else {
            None
        }
    }).collect();
    
    let mut total_force = [0.0, 0.0, 0.0];
    for force_vector in force_vectors {
        total_force[0] += force_vector[0];
        total_force[1] += force_vector[1];
        total_force[2] += force_vector[2];
    }

    let acceleration = [
        total_force[0] / particle.mass,
        total_force[1] / particle.mass,
        total_force[2] / particle.mass,
    ];

    let velocity = [
        particle.velocity[0] + acceleration[0] * TIME_STEP,
        particle.velocity[1] + acceleration[1] * TIME_STEP,
        particle.velocity[2] + acceleration[2] * TIME_STEP,
    ];

    let position = [
        particle.position[0] + particle.velocity[0] * TIME_STEP,
        particle.position[1] + particle.velocity[1] * TIME_STEP,
        particle.position[2] + particle.velocity[2] * TIME_STEP,
    ];

    particle.update_position(position.to_vec());
    particle.update_velocity(velocity.to_vec());
    particle.update_acceleration(acceleration.to_vec());
}

fn force_fluid(particle: &mut Particle, n_bodies: Vec<&Particle>) {}
fn force_electromagnetic(particle: &mut Particle, n_bodies: Vec<&Particle>) {}
// come up with expression for density at a point 
// move particles from high density to low density

fn initialize_particles() -> Vec<Particle> {
    let mut particles: Vec<Particle> = Vec::new();
    for _ in 0..1000 {
        let mass = rand::thread_rng().gen_range(1.0..MASS_EARTH);
        let position = vec![
            rand::thread_rng().gen_range(-100.0..100.0),
            rand::thread_rng().gen_range(-100.0..100.0),
            rand::thread_rng().gen_range(-100.0..100.0),
        ];
        let velocity = vec![
            rand::thread_rng().gen_range(-100.0..100.0),
            rand::thread_rng().gen_range(-100.0..100.0),
            rand::thread_rng().gen_range(-100.0..100.0),
        ];
        let acceleration = vec![
            rand::thread_rng().gen_range(-100.0..100.0),
            rand::thread_rng().gen_range(-100.0..100.0),
            rand::thread_rng().gen_range(-100.0..100.0),
        ];
        particles.push(Particle::new(mass, position, velocity, acceleration));
    }
    return particles;
}

fn three_dimensional_particle_plot(particles: &Vec<Particle>) {
    let root = BitMapBackend::gif(OUT_FILE_NAME, (1024, 768), 100)
        .unwrap()
        .into_drawing_area();
    root.fill(&WHITE).unwrap();

    let min_x = particles.iter().map(|particle| particle.position[0]).min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let max_x = particles.iter().map(|particle| particle.position[0]).max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let min_y = particles.iter().map(|particle| particle.position[1]).min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let max_y = particles.iter().map(|particle| particle.position[1]).max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let min_z = particles.iter().map(|particle| particle.position[2]).min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let max_z = particles.iter().map(|particle| particle.position[2]).max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    
    let mut chart = ChartBuilder::on(&root)
        .caption("3D Scatter Plot", ("sans-serif", 50))
        .build_cartesian_3d(min_x..max_x, min_y..max_y, min_z..max_z)
        .unwrap();

    chart.with_projection(|mut pb| {
        pb.yaw = 0.5;
        pb.scale = 0.9;
        pb.into_matrix()
    });

    chart
        .configure_axes()
        .draw()
        .unwrap();

    chart
        .draw_series(
        PointSeries::of_element(
            particles
                .iter()
                .map(|particle| (particle.position[0], particle.position[1], particle.position[2])),
            5,
            &RED,
            &|c, s, st| {
                return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                    + Circle::new((0, 0), s, st.filled()); // At this point, the new pixel coordinate is established
                    //+ Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
            },

            )
        )
        .unwrap();
}
fn main() {
    let mut particles = initialize_particles();

    for _ in tqdm(0..10) {
        let particles_prior_state = particles.clone();
        for particle in particles.iter_mut() {
            force_gravity(particle, particles_prior_state.iter().collect());
        }
    }
    three_dimensional_particle_plot(&particles);
}
