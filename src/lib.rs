pub type RGB = glm::DVec3;

pub fn zero_dvec3() -> glm::DVec3 {
    glm::dvec3(0., 0., 0.)
}

pub fn mult_add(scalar: f64, a: &glm::DVec3, b: &glm::DVec3) -> glm::DVec3 {
    (glm::to_dvec3(scalar) * (*a)) + (*b)
}

fn point_towards_sphere(
    origin: &glm::DVec3,
    target: &Sphere
) -> glm::DVec3 {
    glm::normalize(
        mult_add(
            -1.,
            origin,
            &target.center
        )
    )
}

#[derive(Debug)]
pub struct Sphere {
    pub center: glm::DVec3,
    pub color: RGB,
    pub radius: f64,
    pub diffuse_coeff: f64,
    pub specular_coeff: f64,
    pub transmission_coeff: f64,
    pub emission_coeff: f64,
    pub refraction_index: f64
}

impl Sphere {
    pub fn new(
        center: glm::DVec3,
        color: RGB,
        radius: f64,
        diffuse_coeff: f64,
        specular_coeff: f64,
        transmission_coeff: f64,
        emission_coeff: f64,
        refraction_index: f64
    ) -> Self {
        Sphere {
            center,
            color,
            radius,
            diffuse_coeff,
            specular_coeff,
            transmission_coeff,
            emission_coeff,
            refraction_index
        }
    }

    #[cfg(test)]
    fn blank_sphere(position: glm::DVec3, radius: f64) -> Self {
        Sphere {
            center: position,
            color: RGB::new(1., 1., 1.),
            radius: radius,
            diffuse_coeff: 0.,
            specular_coeff: 0.,
            transmission_coeff: 0.,
            emission_coeff: 0.,
            refraction_index: 0.
        }
    }
}

#[derive(Debug)]
pub struct Ray {
    origin: glm::DVec3,
    direction: glm::DVec3
}

impl Ray {
    const TOLERANCE: f64 = 1e-7;

    pub fn new(
        origin: glm::DVec3,
        direction: glm::DVec3
    ) -> Self {
        Ray {
            origin: origin,
            direction: direction
        }
    }

    pub fn from_origin(direction: glm::DVec3) -> Self {
        Ray {
            origin: zero_dvec3(),
            direction: direction
        }
    }

    pub fn intersect<'a>(&self, sphere: &'a Sphere) -> Option<(&'a Sphere, f64)> {
        let ray_sphere_comb = mult_add(-1., &self.origin, &sphere.center);
        let ray_sphere_dot = glm::dot(self.direction, ray_sphere_comb);
        let vantage_point = (
            ray_sphere_dot * ray_sphere_dot
        ) - glm::dot(
            ray_sphere_comb,
            ray_sphere_comb
        ) + (
            sphere.radius * sphere.radius
        );
        if vantage_point > 0. {
            let u_sqrt = vantage_point.sqrt();
            let dist = if (ray_sphere_dot - u_sqrt) > Self::TOLERANCE {
                ray_sphere_dot - u_sqrt
            } else {
                ray_sphere_dot + u_sqrt
            };
            if dist >= Self::TOLERANCE {
                Some((sphere, dist))
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn closest_intersection<'a>(&'a self, spheres: &'a Vec<Sphere>) -> Option<Hit> {
        spheres.iter()
            // Filters spheres that the ray intersects
            .filter_map(|sphere| self.intersect(sphere))
            // Gets the sphere with the minimum distance
            .min_by(|(_, a), (_, b)| a.total_cmp(b))
            // Maps sphere to a Hit object
            .map_or(
            None,
            |(sph, dist)| Some(
                Hit::new(
                    sph,
                    self,
                    mult_add(
                        dist,
                        &self.direction,
                        &self.origin
                    )
                )
            )
        )
    }
}

#[derive(Debug)]
pub struct Hit<'a> {
    pub sphere: &'a Sphere,
    ray: &'a Ray,
    intersect_point: glm::DVec3,
    normal: glm::DVec3,
    dot: f64,
    internal_refraction: bool
}

impl<'a> Hit<'a> {
    fn new(sphere: &'a Sphere, ray: &'a Ray, intersect_point: glm::DVec3) -> Self {
        let normal = point_towards_sphere(&intersect_point, sphere);
        let dot = -glm::dot(ray.direction, normal);
        let (normal, dot, internal_refraction) = if dot < 0. {
            (glm::to_dvec3(-1) * normal, -dot, true)
        } else {
            (normal, dot, false)
        };
        Hit {
            sphere,
            ray,
            intersect_point,
            normal,
            dot,
            internal_refraction
        }
    }
}

#[derive(Debug)]
pub struct Scene {
    ambient: glm::DVec3,
    pub spheres: Vec<Sphere>
}

impl Scene {
    pub fn new(
        ambient: RGB,
        spheres: Vec<Sphere>
    ) -> Self {
        Self {
            ambient: ambient,
            spheres: spheres
        }
    }
}

#[derive(Debug)]
pub struct BussinessCardRenderer {
    max_depth: usize,
    image_res: (usize, usize),
    viewing_angle: f64 // in degrees
}

impl BussinessCardRenderer {
    pub fn new(
        max_depth: usize,
        image_res: (usize, usize),
        viewing_angle: f64
    ) -> Self {
        Self {
            max_depth: max_depth,
            image_res: image_res,
            viewing_angle: viewing_angle
        }
    }

    fn ambient(&self, scene: &Scene, hit: &Hit) -> RGB {
        // Ambient color is the color of the ambient light composite with the sphere
        // color and the diffuse coefficient
        glm::to_dvec3(hit.sphere.diffuse_coeff) * hit.sphere.color * scene.ambient
    }

    fn emission(&self, hit: &Hit) -> RGB {
        // Emission is the sphere color times the emission intensity
        glm::to_dvec3(hit.sphere.emission_coeff) * hit.sphere.color
    }

    fn diffuse(&self, scene: &Scene, hit: &Hit) -> RGB {
        // Diffuse is the sum of all the visible lights by the point composite with the
        // sphere color times the diffuse coefficient
        let diffuse_color = scene.spheres.iter()
            .filter_map(
                // For each light
                |light| {
                    // Calculates direction from the point of intersection of initial ray
                    // towards light
                    let vector_towards_light = point_towards_sphere(
                        &hit.intersect_point,
                        light
                    );
                    // Creates new ray pointing towards light
                    let diffuse_ray = Ray::new(
                        hit.intersect_point.clone(),
                        vector_towards_light.clone()
                    );
                    // Calculate by how much intensity is attenuated
                    let attenuation = glm::dot(
                        hit.normal,
                        vector_towards_light
                    );
                    if attenuation > 0. {
                        // Get closest light in path of ray
                        diffuse_ray.closest_intersection(&scene.spheres).and_then(
                            |light_hit| {
                                // If closest intersection is current light and there is still light
                                // remaining after attenuation
                                if std::ptr::eq(light, light_hit.sphere) {
                                    // Get attenuated emission 
                                    Some(glm::to_dvec3(attenuation) * self.emission(&light_hit))
                                } else {
                                    None
                                }
                            }
                        )
                    } else {
                        None
                    }
                }
            )
            // Accumulate intensity of all visible lights
            .fold(
                zero_dvec3(),
                |acum, light_inten| {
                    acum + light_inten
                }
            );
        glm::to_dvec3(hit.sphere.diffuse_coeff) * hit.sphere.color * diffuse_color
    }

    fn specular(&self, scene: &Scene, hit: &Hit, level: usize) -> RGB {
        let specular_origin = &hit.intersect_point;
        let specular_direction = mult_add(
            2. * hit.dot,
            &hit.normal,
            &hit.ray.direction
        );
        glm::to_dvec3(hit.sphere.specular_coeff) * self.trace(
            scene,
            &Ray::new(
                specular_origin.clone(),
                specular_direction
            ),
            level + 1
        )
    }

    fn transmission(&self, scene: &Scene, hit: &Hit, level: usize) -> RGB {
        // Dot has already been flipped by this point
        let material_interface = if hit.internal_refraction {
            1. / hit.sphere.refraction_index
        } else {
            hit.sphere.refraction_index
        };
        let factor = 1. - material_interface * material_interface * (1. - hit.dot * hit.dot);
        if factor > 0. {
            let transmission_origin = &hit.intersect_point;
            let transmission_direction = mult_add(
                material_interface,
                &hit.ray.direction,
                &mult_add(
                    material_interface * hit.dot - factor.sqrt(),
                    &hit.normal,
                    &zero_dvec3()
                )
            );
            glm::to_dvec3(hit.sphere.transmission_coeff) * self.trace(
                scene,
                &Ray::new(
                    transmission_origin.clone(),
                    transmission_direction
                ),
                level + 1
            )
        } else {
            zero_dvec3()
        }
    }

    fn intensity(&self, scene: &Scene, hit: &Hit, level: usize) -> RGB {
        let ambient = self.ambient(scene, hit);
        let emission = self.emission(hit);
        let diffuse = self.diffuse(scene, hit);
        let specular = self.specular(scene, hit, level);
        let transmission = self.transmission(scene, hit, level);
        
        ambient + emission + diffuse + specular + transmission
    }

    pub fn trace<'a>(&'a self, scene: &Scene, ray: &'a Ray, level: usize) -> RGB {
        match level.cmp(&self.max_depth) {
            std::cmp::Ordering::Less => {
                ray.closest_intersection(&scene.spheres).map_or(
                    scene.ambient,
                    |hit| self.intensity(&scene, &hit, level)
                )
            },
            _ => zero_dvec3()
        }
    }

    pub fn render_ppm<T>(&self, io: &mut T, scene: &Scene) where T: std::io::Write {
        let (w, h) = self.image_res;
        let (w, h) = (w as isize, h as isize);

        writeln!(io, "P3").unwrap();
        writeln!(io, "{} {}", w, h).unwrap();
        writeln!(io, "255").unwrap();
        let magnetude = (w as f64) / 2. / glm::radians(self.viewing_angle / 2.).tan();
        for z in (0..h).rev() {
            for x in 0..w {
                let pixel_ray = Ray::from_origin(
                    glm::normalize(
                        glm::dvec3(
                            (x - w/2) as f64,
                            magnetude,
                            (1 + z - h/2) as f64,
                        )
                    )
                );
                let color = mult_add(
                    255.,
                    &self.trace(scene, &pixel_ray, 0),
                    &zero_dvec3()
                );
                writeln!(
                    io,
                    "{} {} {}",
                    color.x.clamp(0., 255.) as u8,
                    color.y.clamp(0., 255.) as u8,
                    color.z.clamp(0., 255.) as u8
                ).unwrap();
            }
        }
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use crate::*;

    #[test]
    fn intersect_direct_hit() {
        let sphere = Sphere::blank_sphere(glm::dvec3(0., 0., 5.), 1.);
        let ray = Ray::from_origin(glm::dvec3(0., 0., 1.));
        let hit = ray.intersect(&sphere).unwrap();

        assert_eq!(hit.1, 4.);
    }

    #[test]
    fn intersect_direct_hit_far() {
        let sphere = Sphere::blank_sphere(glm::dvec3(0., 0., 50000.), 1.);
        let ray = Ray::from_origin(glm::dvec3(0., 0., 1.));
        let hit = ray.intersect(&sphere).unwrap();

        assert_eq!(hit.1, 49999.);
    }

    #[test]
    #[should_panic]
    fn intersect_glancing_hit() {
        let sphere = Sphere::blank_sphere(glm::dvec3(0., 1., 5.), 1.);
        let ray = Ray::from_origin(glm::dvec3(0., 0., 1.));
        ray.intersect(&sphere).unwrap();
    }

    #[test]
    #[should_panic]
    fn intersect_miss() {
        let sphere = Sphere::blank_sphere(glm::dvec3(1., 1., 5.), 1.);
        let ray = Ray::from_origin(glm::dvec3(0., 0., 1.));
        ray.intersect(&sphere).unwrap();
    }

    #[test]
    fn closest_intersection() {
        let spheres = vec![
            Sphere::blank_sphere(glm::dvec3(0., 0., 5.), 1.),
            Sphere::blank_sphere(glm::dvec3(0., 0., 4.), 1.),
        ];
        let ray = Ray::from_origin(glm::dvec3(0., 0., 1.));
        let hit = ray.closest_intersection(&spheres).unwrap();

        assert!(std::ptr::eq(hit.sphere, &spheres[1]));
        assert_eq!(hit.intersect_point, glm::dvec3(0., 0., 3.));
    }

    #[test]
    fn closest_intersection2() {
        let spheres = vec![
            Sphere::blank_sphere(glm::dvec3(0., 0., 4.), 1.),
            Sphere::blank_sphere(glm::dvec3(0., 0., 5.), 1.),
        ];
        let ray = Ray::from_origin(glm::dvec3(0., 0., 1.));
        let hit = ray.closest_intersection(&spheres).unwrap();

        assert!(std::ptr::eq(hit.sphere, &spheres[0]));
        assert_eq!(hit.intersect_point, glm::dvec3(0., 0., 3.));
    }

    #[test]
    fn closest_intersection3() {
        let mut spheres: Vec<Sphere> = vec![];
        for i in 5..20 {
            spheres.push(
                Sphere::blank_sphere(glm::dvec3(0., 0., i as f64), 1.)
            );
        }
        spheres.push(Sphere::blank_sphere(glm::dvec3(0., 0., 4.), 1.));
        let ray = Ray::from_origin(glm::dvec3(0., 0., 1.));
        let hit = ray.closest_intersection(&spheres).unwrap();

        assert!(std::ptr::eq(hit.sphere, spheres.last().unwrap()));
        assert_eq!(hit.intersect_point, glm::dvec3(0., 0., 3.));
    }

    #[test]
    fn ray_directions() {
        use glm::is_close_to;
        const TOLERANCE: f64 = f64::EPSILON;

        const SIZE: usize = 700;
        const VIEWING_ANGLE: f64 = 25.;

        let (w, h) = (SIZE, SIZE);
        let (wi, hi) = (w as isize, h as isize);
        let wf = w as f64;

        let magnetude = (wf) / 2. / glm::radians(VIEWING_ANGLE / 2.).tan();
        let c_iter = 0..(wi * hi);
        let rust_iter = ((0..hi).rev()).cartesian_product(0..wi);
        for (yx, (z, x)) in c_iter.zip(rust_iter) {
            let c_ray_direction = glm::dvec3(
                (yx % wi - wi / 2) as f64, 
                magnetude,
                (wi / 2 - yx / wi) as f64
            );
            let rust_ray_direction = glm::dvec3(
                (x - wi/2) as f64,
                magnetude,
                (1 + z - hi/2) as f64,
            );

            glm::assert_close_to!(c_ray_direction, rust_ray_direction, TOLERANCE);
        }
    }
}