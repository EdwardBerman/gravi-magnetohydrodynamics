const GRAVITATIONAL_CONSTANT: f64 = 6.67408e-11;

struct Particle {
    mass: f64,
    position: Vec<f64>,
}

impl Particle {
    fn new(mass: f64, position: Vec<f64>) -> Particle {
        Particle {
            mass,
            position,
        }
    }
}

impl Particle {
    fn get_closest_neighbors(&self, neighbors: Vec<Particle>) -> Vec<Particle> {
        let mut closest_neighbors: Vec<Particle> = Vec::new();
        for neighbor in neighbors {
            let distance = (self.position[0] - neighbor.position[0]).powf(2.0);
            if distance < 1.0 {
                closest_neighbors.push(neighbor);
            }
        }
        return closest_neighbors;
    }
}

impl Particle {
    fn update_position(&mut self, position: Vec<f64>) {
        self.position = position;
    }
}

fn force_gravity(particle: &mut Particle, n_bodies: Vec<Particle>) -> &mut Particle {
    for neighbor in particle.get_closest_neighbors(n_bodies) {
        let distance = (particle.position[0] - neighbor.position[0]).powf(2.0);
        let force = GRAVITATIONAL_CONSTANT * (particle.mass * neighbor.mass) / distance;
        let acceleration = force / particle.mass;
        let velocity = acceleration * 1.0;
        let position = velocity * 1.0;
        // vectorize force, velocity, and position
        let force = vec![force];
        let velocity = vec![velocity];
        let position = vec![position];
        particle.update_position(position);
    }
    return particle;
}

fn main() {
    println!("Hello, world!");
}
