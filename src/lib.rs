//This code has not yet been tested. It may not be correct.

use std::collections::HashMap;

use nalgebra::{Matrix3, Matrix4, Vector3, Unit};
use noise::{NoiseFn, Simplex};
use num_complex::Complex64;
use rand::prelude::*;

type PoincarePoint = Vector3<f64>;
type PoincareHalfPlanePoint = Complex64;
type MinkowskiPoint = Vector3<f64>;
type BarycentricPoint = Vector3<f64>;

pub struct H2xE {
    noise: Simplex,
    rng: ThreadRng,
}

impl H2xE {
    pub fn new() -> Self {
        Self {
            noise: Simplex::new(0), //TODO: Seed mechanism
            rng: rand::thread_rng(),
        }
    }

    pub fn hyperboloid_matrix() -> Matrix4<f64> {
        Matrix4::new(
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, -1.0,
        )
    }

    /// Convert a point from Poincare disk model to Minkowski hyperboloid model.
    /*pub fn poincare_to_hyperboloid(&self, poincare_point: &Vector3<f64>) -> Vector3<f64> {
        let r2 = poincare_point.y.powi(2) + poincare_point.z.powi(2);
        let x = (1.0 + r2) / (1.0 - r2);
        let y = 2.0 * poincare_point.y / (1.0 - r2);
        let z = 2.0 * poincare_point.z / (1.0 - r2);
        Vector3::new(x, y, z)
    }*/
    pub fn poincare_to_hyperboloid(&self, point: &PoincarePoint) -> MinkowskiPoint {
        let r2 = point.x * point.x + point.y * point.y + point.z * point.z;
        let x = (r2 + 1.0) / (1.0 - r2);
        let y = 2.0 * point.x / (1.0 - r2);
        let z = 2.0 * point.y / (1.0 - r2);
        Vector3::new(x, y, z)
    }

    /// Convert a point from Minkowski hyperboloid model to Poincare disk model.
    /*pub fn hyperboloid_to_poincare(&self, hyperboloid_point: &Vector3<f64>) -> Vector3<f64> {
        let x = hyperboloid_point.x;
        let y = hyperboloid_point.y / (x + 1.0);
        let z = hyperboloid_point.z / (x + 1.0);
        Vector3::new(x, y, z)
    }*/
    pub fn hyperboloid_to_poincare(&self, point: &MinkowskiPoint) -> PoincarePoint {
        let x = point.y / (point.x + 1.0);
        let y = point.z / (point.x + 1.0);
        let z = point.z;
        Vector3::new(x, y, z)
    }

    /// Calculate the distance between two points in the Minkowski hyperboloid model.
    /*pub fn distance(&self, point1: &Vector3<f64>, point2: &Vector3<f64>) -> f64 {
        let inner_product = point1.x * point2.x - point1.y * point2.y - point1.z * point2.z;
        let cosh_dist = inner_product / ((point1.x - 1.0) * (point2.x - 1.0));

        cosh_dist.acosh()
    }*/
    pub fn distance(&self, point_a: &MinkowskiPoint, point_b: &MinkowskiPoint) -> f64 {
        let delta = point_a - point_b;
        let delta_norm = delta.norm();
        let angle = ((point_a.x * point_b.x) - (point_a.y * point_b.y) - (point_a.z * point_b.z)).acos();
        angle.sinh() * delta_norm
    }

    pub fn direction(&self, point_a: &MinkowskiPoint, point_b: &MinkowskiPoint) -> Vector3<f64> {
        let tangent = point_b - point_a;
        let weight = (point_a.x + 1.0) / point_a.x;
        tangent / weight
    }

    /// Calculate the geodesic between two points in the Minkowski hyperboloid model.
    /*pub fn geodesic(&self, point1: &MinkowskiPoint, point2: &MinkowskiPoint, t: f64) -> MinkowskiPoint {
        let dist = self.distance(point1, point2);
        let direction = (point2 - point1).normalize();

        // Calculate the exponential map of the direction vector multiplied by the distance and t (0 <= t <= 1)
        let mut exp_map = Matrix3::identity();
        let factor = dist * t;
        exp_map[(1, 1)] = factor.cos();
        exp_map[(2, 2)] = factor.cos();
        exp_map[(1, 2)] = factor.sin();
        exp_map[(2, 1)] = -factor.sin();

        exp_map * point1
    }*/
    pub fn geodesic(&self, point1: &MinkowskiPoint, point2: &MinkowskiPoint, t: f64) -> MinkowskiPoint {
        let dist = self.distance(point1, point2);
        let direction = (point2 - point1).normalize();
    
        // Calculate the exponential map of the direction vector multiplied by the distance and t (0 <= t <= 1)
        let factor = dist * t;
        let cosh_factor = factor.cosh();
        let sinh_factor = factor.sinh();
    
        // Update the Minkowski point using the exponential map
        MinkowskiPoint::new(point1.x, point1.y * cosh_factor + direction.y * sinh_factor, point1.z * cosh_factor + direction.z * sinh_factor)
    }
    

    /// Update the position of an object in the Minkowski hyperboloid model based on its velocity and time step.
    pub fn update_position(
        &self,
        position: &MinkowskiPoint,
        velocity: &Vector3<f64>,
        time_step: f64,
    ) -> MinkowskiPoint {
        // Calculate the exponential map of the velocity vector in the XY plane (H2 part) multiplied by the time step
        let mut exp_velocity_h2 = Matrix3::identity();
        let factor_h2 = (velocity.y.powi(2) + velocity.z.powi(2)).sqrt() * time_step;
        exp_velocity_h2[(1, 1)] = factor_h2.cos();
        exp_velocity_h2[(2, 2)] = factor_h2.cos();
        exp_velocity_h2[(1, 2)] = factor_h2.sin();
        exp_velocity_h2[(2, 1)] = -factor_h2.sin();

        // Update the position in the XY plane (H2 part) by applying the exponential map to the current position
        let new_position_h2 = exp_velocity_h2 * position;

        // Update the position in the Z-axis (E part) linearly
        let new_position_z = position.z + velocity.z * time_step;

        // Combine the updated H2 and E parts to get the new position in H2xE space
        Vector3::new(new_position_h2.x, new_position_h2.y, new_position_z)
   
    }

    /// Perform a rotation in H2xE space.
    /// `point` is the Minkowski hyperboloid point to be rotated.
    /// `angle_h2` is the angle of rotation in the H2 plane (hyperbolic plane) in radians.
    /// `angle_e` is the angle of rotation in the E component (Euclidean dimension) in radians.
    /*pub fn rotate(&self, point: &Vector3<f64>, angle_h2: f64, angle_e: f64) -> Vector3<f64> {
        let rotation_h2 = Matrix3::new(
            1.0, 0.0, 0.0,
            0.0, angle_h2.cos(), angle_h2.sin(),
            0.0, -angle_h2.sin(), angle_h2.cos()
        );

        let rotation_e = Matrix3::new(
            angle_e.cos(), -angle_e.sin(), 0.0,
            angle_e.sin(), angle_e.cos(), 0.0,
            0.0, 0.0, 1.0
        );

        let rotated_point = rotation_h2 * point;
        rotation_e * rotated_point
    }*/
    pub fn rotate(&self, point: &MinkowskiPoint, angle_h2: f64, angle_e: f64) -> MinkowskiPoint {
        let poincare_point = self.hyperboloid_to_poincare(point);
    
        // Perform 2D Euclidean rotation for the H2 part (Poincare disk model)
        let rotated_x = poincare_point.x * angle_h2.cos() - poincare_point.y * angle_h2.sin();
        let rotated_y = poincare_point.x * angle_h2.sin() + poincare_point.y * angle_h2.cos();
    
        // Perform 1D rotation for the E part (Euclidean dimension)
        let rotated_z = poincare_point.z * angle_e.cos();
    
        let rotated_poincare_point = Vector3::new(rotated_x, rotated_y, rotated_z);
        self.poincare_to_hyperboloid(&rotated_poincare_point)
    }
    

    /// Sample the 2D simplex noise function at a given point in the Minkowski hyperboloid model.
    pub fn sample_noise_at(&self, point: &MinkowskiPoint) -> f64 {
        let poincare_point = self.hyperboloid_to_poincare(point);
        self.noise.get([poincare_point.x, poincare_point.y])
    }

    /// Generate random barycentric coordinates within a cell of the hyperbolic tessellation.
    pub fn random_barycentric_coordinates(&mut self) -> BarycentricPoint {
        let u = self.rng.gen::<f64>();
        let v = self.rng.gen::<f64>();

        if u + v > 1.0 {
            Vector3::new(1.0 - u, 1.0 - v, 1.0 - (1.0 - u) - (1.0 - v))
        } else {
            Vector3::new(u, v, 1.0 - u - v)
        }
    }

    /// Convert barycentric coordinates within a cell of the hyperbolic tessellation to Cartesian coordinates in the Poincare disk model.
    pub fn barycentric_to_poincare(&self, barycentric: BarycentricPoint, cell_vertices: [Vector3<f64>; 3]) -> PoincarePoint {
        let (u, v, w) = (barycentric.x, barycentric.y, barycentric.z);
        let poincare_vertices = [
            self.hyperboloid_to_poincare(&cell_vertices[0]),
            self.hyperboloid_to_poincare(&cell_vertices[1]),
            self.hyperboloid_to_poincare(&cell_vertices[2]),
        ];

        u * poincare_vertices[0] + v * poincare_vertices[1] + w * poincare_vertices[2]
    }
}

pub trait PoincarePointMethods {
    fn to_upper_half_plane(&self) -> PoincareHalfPlanePoint;
    fn from_upper_half_plane(z: PoincareHalfPlanePoint) -> Self;
}

impl PoincarePointMethods for PoincarePoint {
    fn to_upper_half_plane(&self) -> PoincareHalfPlanePoint {
        let x = self.x;
        let y = self.y;
        let denominator = x * x + (1.0 - y) * (1.0 - y);
        let real = (2.0 * x) / denominator;
        let imaginary = (1.0 - x * x - y * y) / denominator;

        Complex64::new(real, imaginary)
    }

    fn from_upper_half_plane(z: PoincareHalfPlanePoint) -> Self {
        let x = (4.0 * z.re) / ((z.re * z.re) + (z.im * z.im));
        let y = (z.re * z.re + z.im * z.im - 1.0) / ((z.re * z.re) + (z.im * z.im));

        Vector3::new(x, y, 0.0)
    }
}

//#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ChunkId(Vec<u8>);

pub struct Chunk {
    pub id: ChunkId,
    //pub data: Vec<Vec<f32>>, //Just an example
}

pub struct ChunkManager {
    pub chunks: HashMap<ChunkId, Chunk>,
    //pub chunk_size: f64,
    pub max_depth: usize,
}

impl ChunkManager {
    //pub fn new(chunk_size: f32, max_depth: usize) -> Self {
    pub fn new(max_depth: usize) -> Self {
        Self {
            chunks: HashMap::new(),
            //chunk_size,
            max_depth,
        }
    }

    pub fn get_chunk(&mut self, id: ChunkId) -> &Chunk {
        if !self.chunks.contains_key(&id) {
            // Generate chunk data here
            //let data = vec![vec![0.0; self.chunk_size as usize]; self.chunk_size as usize];
            //let chunk = Chunk { id, data };
            let chunk = Chunk { id: id.clone() };
            self.chunks.insert(id.clone(), chunk);
        }
        self.chunks.get(&id).unwrap()
    }

    pub fn get_chunk_for_position(&mut self, position: &PoincarePoint) -> &Chunk {
        //let chunk_id = self.to_chunk_id(position, self.chunk_size);
        let chunk_id = self.to_chunk_id(position); //TODO: Proper max_depth implementation
        self.get_chunk(chunk_id)
    }

    fn to_chunk_id(&self, poincare_point: &PoincarePoint) -> ChunkId {
        let point = poincare_point.to_upper_half_plane();
        let mut chunk_id = Vec::with_capacity(self.max_depth);
        let mut current_point = point;
        let mut depth = 0;
    
        let f1 = |z: PoincareHalfPlanePoint| (z - 1.0) / z;
        let f2 = |z: PoincareHalfPlanePoint| z + 2.0;
        let f3 = |z: PoincareHalfPlanePoint| z + 1.0;
    
        while depth < self.max_depth {
            let a = f1(current_point).norm_sqr();
            let b = f2(current_point).norm_sqr();
            let c = f3(current_point).norm_sqr();
    
            let min_norm = a.min(b).min(c);
    
            if a == min_norm {
                chunk_id.push(0);
                current_point = f1(current_point);
            } else if b == min_norm {
                chunk_id.push(1);
                current_point = f2(current_point);
            } else {
                chunk_id.push(2);
                current_point = f3(current_point);
            }
    
            depth += 1;
        }
    
        ChunkId(chunk_id)
    }

    fn from_chunk_id(&self, chunk_id: ChunkId) -> PoincarePoint {
        let f1_inv = |z: PoincareHalfPlanePoint| z / (z - 1.0);
        let f2_inv = |z: PoincareHalfPlanePoint| z - 2.0;
        let f3_inv = |z: PoincareHalfPlanePoint| z - 1.0;
    
        let mut current_point = PoincareHalfPlanePoint::new(0.0, 1.0);
    
        for &operation in chunk_id.0.iter().take(self.max_depth) {
            match operation {
                0 => current_point = f1_inv(current_point),
                1 => current_point = f2_inv(current_point),
                2 => current_point = f3_inv(current_point),
                _ => panic!("Invalid operation in chunk_id"),
            }
        }
    
        PoincarePoint::from_upper_half_plane(current_point)
    }
    
}