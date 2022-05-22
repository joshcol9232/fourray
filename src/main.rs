use nalgebra as na;

use std::io::Write;
use na::{Point4, Vector4, Vector3};
use image::{DynamicImage, GenericImage, Rgba};

const DIMS: [u32; 2] = [800, 800];

pub type Colour = Vector3<f64>;


struct Ray {
    orig: Point4<f64>,
    dir: Vector4<f64>,
    pixel: [u32; 2],
    colour: Colour,
}

impl Ray {
    pub fn new(o: Point4<f64>, d: Vector4<f64>, i: u32, j: u32) -> Self {
        Self {
            orig: o,
            dir: d,
            pixel: [i, j],
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

/* 
fn cross4D(a: &Vector4<f64>, b: &Vector4<f64>, c: &Vector4<f64>) -> Vector4<f64> {
    // From: Raytracing 4D fractals, visualising the four dimensional
    // properties of the Julia set; Torkel Odegaard, Joakim Wennergren
    Vector4::new(
        a[1] * (b[3] * c[2] - b[2] * c[3]) + b[1] * (a[2] * c[3] - a[3] * c[2]) + c[1] * (a[3] * b[2] - a[2] * b[3]),
        a[0] * (b[2] * c[3] - b[3] * c[2]) + b[0] * (a[3] * c[2] - a[2] * c[3]) + c[0] * (a[2] * b[3] - a[3] * b[2]),
        a[0] * (b[3] * c[1] - b[1] * c[3]) + b[0] * (a[1] * c[3] - a[3] * c[1]) + c[0] * (a[3] * b[1] - a[1] * b[3]),
        a[0] * (b[1] * c[2] - b[2] * c[1]) + b[0] * (a[2] * c[1] - a[1] * c[2]) + c[0] * (a[1] * b[2] - a[2] * b[1]),
    )
}
*/

fn get_rgba(c: Colour) -> Rgba<u8> {
    Rgba([
        (c[0] * 255.0).ceil() as u8,
        (c[1] * 255.0).ceil() as u8,
        (c[2] * 255.0).ceil() as u8,
        255,
    ])
}

fn spawn_rays(origin: Point4<f64>, viewdir: Vector4<f64>, upvec: Vector4<f64>, limbovec: Vector4<f64>, fov: f64, xmax: u32, ymax: u32) -> Vec<Ray> {
    /* Idea: Viewport is plane of dims 1 x 1
    Scan along viewport in dx = 1/X, dy = 1/Y, starting at x = -1/2, y = -1/2
    Offset viewport from origin by some distance `d`.
    `d` is then d = 1/2tan(fov/2).
    For each pixel on the viewport, find "surface normal" from origin as if a sphere.
    Spawn rays with origin on the plane, with direction from sphere surface of dist to origin
    This is just normalised origin? not sure

    Source of up vector & limbo vector idea: Raytracing 4D fractals, visualising the four dimensional
    properties of the Julia set; Torkel Odegaard, Joakim Wennergren

    NOTE: For now, just debugging with known directions
    */
    
    let viewport_dist = 1.0/(2.0 * (fov/2.0).tan()); // This is a point on our plane at (x=0, y=0)
    let mut x = -0.5;  // In plane of viewport
    let mut y = -0.5;

    let dx = 1.0/xmax as f64;
    let dy = dx;

    let mut rays = Vec::with_capacity(xmax as usize * ymax as usize);

    let plane_centre = origin + viewport_dist * Vector4::new(0.0, 0.0, viewport_dist, 0.0);  // View along +ve z axis
    
    let mut j = 0;

    while j < ymax { // scanlines
        let mut i = 0;
        x = -0.5;
        while i < xmax {
            let mut o = plane_centre.clone();
            o[0] += x;
            o[1] += y;

            let mut d = (o - origin);
            d /= d.norm();
            let r = Ray::new(o, d, i, j);

            rays.push(r);
            // println!("{} {}", i, j);

            x += dx;
            i += 1;
        }

        y += dy;
        j += 1;
    }

    rays
}

fn main() {
    let mut img = DynamicImage::new_rgb8(DIMS[0], DIMS[1]);

    let rays = spawn_rays(
        Point4::new(0.0, 0.0, 0.0, 0.0),
        Vector4::new(0.0, 0.0, 1.0, 0.0),  // NOTE: Doesn't make any difference atm
        Vector4::new(0.0, 1.0, 0.0, 0.0),
        Vector4::new(0.0, 0.0, 0.0, 1.0),  // Effectively makes it 3D
        90.0,
        DIMS[0],
        DIMS[1],
    );

    let objects = vec![
        Sphere::new(Point4::new(0.0, 0.0, 1.0, 0.0), 1.0, Colour::new(0.0, 1.0, 0.0)),
    ];

    for (i, ray) in rays.iter().enumerate() {
        // println!("Marching ray {}/{}.", i, rays.len());
        let hitcolour: Colour = match ray.march(0.01, &objects, 1.0, 10.0) {
            Some(hit) => hit.colour,
            None => Colour::new(0.0, 0.0, 0.0),
        };

        img.put_pixel(ray.pixel[0], ray.pixel[1], get_rgba(hitcolour));
    }


    img.save("out.png").unwrap();
}