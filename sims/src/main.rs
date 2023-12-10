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
    let mut total_force = vec![0.0, 0.0, 0.0];
    for neighbor in particle.get_closest_neighbors(n_bodies) {
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

            total_force[0] += force_vector[0];
            total_force[1] += force_vector[1];
            total_force[2] += force_vector[2];
        }
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

fn main() {
    let mut particle_a = Particle::new(
        500.0,
        vec![0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0],
    );
    let mut particle_b = Particle::new(
        500.0,
        vec![1.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0],
    );

    // same loop as above but only print particle a every 10 iterations
    for _ in 0..500 {
        force_gravity(&mut particle_a, vec![&particle_b]);
        force_gravity(&mut particle_b, vec![&particle_a]);
    }
}
