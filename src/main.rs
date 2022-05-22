use nalgebra as na;

use std::io::Write;
use na::{Point4, Vector4, Vector3};
use image::DynamicImage;

const DIMS: [u32; 2] = [800, 600];

pub type Colour = Vector3<f64>;


struct Ray {
    orig: Point4<f64>,
    dir: Vector4<f64>,
    colour: Colour,
}

impl Ray {
    pub fn new(o: Point4<f64>, d: Vector4<f64>) -> Self {
        Self {
            orig: o,
            dir: d,
            colour: Colour::new(1.0, 1.0, 1.0),
        }
    }

    pub fn march(&self, dt: f64, objects: &[Sphere], nearplane: f64, farplane: f64) -> Option<RayHit> {
        // March ray forwards in increments of dt, from nearplane up to farplane
        let mut t: f64 = nearplane;
        let mut raypos: Point4<f64> = self.orig.clone();

        // List of hits found (all objects must be checked, the nearest in `t` is correct hit).
        let mut current_hits = vec![0.0; objects.len()];

        while t < farplane {
            for (i, obj) in objects.iter().enumerate() {

                if obj.is_point_inside(&raypos) {
                    return Some(RayHit {        // Break as an object has been hit
                        pos: raypos.clone(),
                        colour: self.colour.component_mul(&obj.colour),
                    })
                }
            }
            // Ray defined by P(t) = origin + t * dir
            raypos += self.dir * dt;
            t += dt;
        }

        None
    }
}

struct RayHit {
    pos: Point4<f64>,
    colour: Colour,
}

struct Sphere {
    pos: Point4<f64>,
    rad: f64,
    colour: Colour,
}

impl Sphere {
    pub fn new(pos: Point4<f64>, rad: f64, colour: Colour) -> Self {
        Self {
            pos,
            rad,
            colour,
        }
    }

    pub fn is_point_inside(&self, p: &Point4<f64>) -> bool {
        (p - self.pos).norm_squared() < self.rad.powi(2)
    }

    pub fn surface_normal(&self, p: &Point4<f64>) -> Vector4<f64> {
        let mut v = p - self.pos;
        v /= v.norm();
        v
    }
}


fn spawn_rays(viewdir: Vector4<f64>, upvec: Vector4<f64>, fov: f64, xmax: u32, ymax: u32) -> Vec<Ray> {
    /* Idea: Viewport is plane of dims 1 x 1
    Scan along viewport in dx = 1/X, dy = 1/Y, starting at x = -1/2, y = -1/2
    Offset viewport from origin by some distance `d`.
    `d` is then d = 1/2tan(fov/2).
    For each pixel on the viewport, find "surface normal" from origin as if a sphere.
    Spawn rays with origin on the plane, with direction from sphere surface of dist to origin
    This is just normalised origin? not sure

    Source of up vector idea: Raytracing 4D fractals, visualising the four dimensional
    properties of the Julia set; Torkel Odegaard, Joakim Wennergren
    */
    
    let viewport_dist = 1.0/(2.0 * (fov/2.0).tan()); // This is a point on our plane at (x=0, y=0)
    let mut x = -0.5;  // In plane of viewport
    let mut y = -0.5;

    let dx = 1.0/xmax;
    let dy = 1.0/ymax;
    
    while y < 0.5 { // scanlines
        while x < 0.5 {


            x += dx;
        }

        y += dy;
    }
}

fn main() {
    use std::time::Instant;

    let mut stdout = std::io::stdout();

    let img = DynamicImage::new_rgb8(DIMS[0], DIMS[1]);

    for j in 0..DIMS[1] {
        print!("\rScanlines remaining: {}/{}.", j, DIMS[1]);
        stdout.flush().unwrap();
        for i in 0..DIMS[0] {
            
        }
    }


    img.save("out.png").unwrap();
}