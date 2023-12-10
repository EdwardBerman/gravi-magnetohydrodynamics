use rand::Rng;
use rayon::prelude::*;

const GRAVITATIONAL_CONSTANT: f64 = 6.67408e-11; // m^3 kg^-1 s^-2
const TIME_STEP: f64 = 1.0; // s (for now)

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

fn initialize_particles() -> Vec<Particle> {
    let mut particles: Vec<Particle> = Vec::new();
    for _ in 0..100 {
        let mass = rand::thread_rng().gen_range(1.0..100.0);
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
fn main() {
    let mut particles = initialize_particles();
    let particles_original_state = particles.clone();
    
    for _ in 0..100 {
        for particle in particles.iter_mut() {
            force_gravity(particle, particles_original_state.iter().collect());
        }
    }
}
