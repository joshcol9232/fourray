use nalgebra as na;

use std::{io::Write, mem::Discriminant};
use std::boxed::Box;
use na::{Vector4, Vector3, Normed};
use image::{DynamicImage, GenericImage, Rgba};

const DIMS: [u32; 2] = [1200, 1200];

pub type Colour = Vector3<f64>;


struct RayHit {
    pos: Vector4<f64>,
    colour: Colour,
    shade: f64,
}


struct Ray {
    orig: Vector4<f64>,
    dir: Vector4<f64>,
    pixel: [u32; 2],
    colour: Colour,
}

impl Ray {
    pub fn new(o: Vector4<f64>, d: Vector4<f64>, i: u32, j: u32) -> Self {
        Self {
            orig: o,
            dir: d,
            pixel: [i, j],
            colour: Colour::new(1.0, 1.0, 1.0),
        }
    }

    #[inline]
    fn at(&self, t: f64) -> Vector4<f64> {      // Ray equation
        self.orig + t * self.dir
    }

    pub fn trace(&self, objects: &[Box<dyn Object>]) -> Option<RayHit> {
        for (i, obj) in objects.iter().enumerate() {
            let intersect = obj.ray_intersection(self); // Returns depth at intersection
            // println!("Intersect: {}", intersect);
            if intersect > -1.0 {  // Has hit
                let raypos = self.at(intersect);

                // APPLY BASIC SHADING
                // TODO: Make shading better. Just apply some shading from far away on +ve all axis for now
                let mut norm_line_to_light = Vector4::new(100.0, -100.0, -50.0, 0.0) - raypos;
                norm_line_to_light /= norm_line_to_light.norm();


                let shade = obj
                    .surface_normal(&raypos)
                    .dot(&norm_line_to_light)
                    .max(0.0);
                    
                // println!("Shade: {}", shade);

                return Some(RayHit {        // Break as an object has been hit
                    pos: raypos.clone(),
                    colour: self.colour.component_mul(obj.colour()),
                    shade,
                })
            }
        }

        None
    }
}


trait Object {
    fn ray_intersection(&self, ray: &Ray) -> f64;
    fn colour(&self) -> &Colour;
    fn surface_normal(&self, p: &Vector4<f64>) -> Vector4<f64>;
}

struct Sphere {
    pos: Vector4<f64>,
    rad: f64,
    colour: Colour,
}

impl Sphere {
    pub fn new(pos: Vector4<f64>, rad: f64, colour: Colour) -> Self {
        Self {
            pos,
            rad,
            colour,
        }
    }
}
impl Object for Sphere {
    fn ray_intersection(&self, ray: &Ray) -> f64 {
        // https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
        let diff = ray.orig - self.pos;
        let b = 2.0 * ray.dir.dot(&diff);
        let c = diff.norm_squared() - self.rad*self.rad;
        // a = 1 because a is the dot product of the ray direction with itself, so 1 squared.

        let discriminant = b.powi(2) - 4.0 * c;
        if discriminant >= 0.0 {
            // Now find location of hit
            (-b - discriminant.sqrt() ) / 2.0
        } else {
            -1.0
        }
    }

    fn surface_normal(&self, p: &Vector4<f64>) -> Vector4<f64> {
        let mut v = p - self.pos;
        v /= v.norm();
        v
    }

    fn colour(&self) -> &Colour {
        &self.colour
    }
}

/*
struct Floor {
    y: f64,
    colour: Colour,
}

impl Object for Floor {
    fn is_point_inside(&self, p: &Vector4<f64>) -> bool {
        p[1] > self.y       // Below floor (+ve y is down)
    }

    fn surface_normal(&self, p: &Vector4<f64>) -> Vector4<f64> {
        Vector4::new(0.0, 1.0, 0.0, 0.0)
    }

    fn colour(&self) -> &Colour {
        &self.colour
    }
}
*/

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

fn spawn_rays(origin: Vector4<f64>, viewdir: Vector4<f64>, upvec: Vector4<f64>, limbovec: Vector4<f64>, fov: f64, xmax: u32, ymax: u32) -> Vec<Ray> {
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
    
    let viewport_dist = 1.0/(2.0 * (fov/4.0).tan()); // This is a point on our plane at (x=0, y=0)
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

            let mut d = o - origin;
            d /= d.norm();
            let r = Ray::new(o, d, i, j);

            rays.push(r);

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
        Vector4::new(0.0, 0.0, 0.0, 0.0),
        Vector4::new(0.0, 0.0, 1.0, 0.0),  // NOTE: Doesn't make any difference atm
        Vector4::new(0.0, 1.0, 0.0, 0.0),
        Vector4::new(0.0, 0.0, 0.0, 1.0),  // Effectively makes it 3D
        90.0,
        DIMS[0],
        DIMS[1],
    );

    let objects: Vec<Box<dyn Object>> = vec![
        Box::new(Sphere::new(Vector4::new(1.2, 0.3, 3.0, 0.0), 0.5, Colour::new(0.0, 1.0, 0.0))),
        Box::new(Sphere::new(Vector4::new(-1.2, 0.3, 3.0, 0.0), 0.5, Colour::new(1.0, 0.0, 0.0))),
        // Box::new(Floor { y: 1.5, colour: Colour::new(1.0, 0.0, 1.0) }),
    ];


    for (i, ray) in rays.iter().enumerate() {
        if i as u32 % DIMS[1] == 0 {
            println!("Done {}/{} scanlines.", i as u32/DIMS[1], DIMS[1]);
        }
        
        let hitcolour: Colour = match ray.trace(&objects) {
            Some(hit) => hit.colour * hit.shade,
            None => Colour::new(0.8, 0.8, 0.8),
        };

        img.put_pixel(ray.pixel[0], ray.pixel[1], get_rgba(hitcolour));
    }


    img.save("out.png").unwrap();
}