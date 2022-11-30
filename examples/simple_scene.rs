use heckbert_renderer::*;

fn main() {
    const BOUNCES: usize = 5;
    const SIZE: usize = 700;
    const VIEWING_ANGLE: f64 = 25.;

    let renderer = BussinessCardRenderer::new(
        BOUNCES,
        (SIZE, SIZE),
        VIEWING_ANGLE
    );
    let scene = Scene::new(
        RGB::new(0.02, 0.02, 0.02),
        vec![
            Sphere::new(
                glm::dvec3(0., 6., 0.5),
                RGB::new(1., 1., 1.),
                0.9,
                0.05, 0.2, 0.85, 0.,
                1.7
            ),
            Sphere::new(
                glm::dvec3(-1., 8., -0.5),
                RGB::new(1., 0.5, 0.2),
                1.,
                0.7, 0.3, 0., 0.05,
                1.2
            ),
            Sphere::new(
                glm::dvec3(1., 8., -0.5),
                RGB::new(0.1, 0.8, 0.8),
                1.,
                0.3, 0.7, 0., 0.,
                1.2
            ),
            Sphere::new(
                glm::dvec3(3., -6., 15.),
                RGB::new(1., 0.8, 1.),
                7.,
                0., 0., 0., 0.6,
                1.5
            ),
            Sphere::new(
                glm::dvec3(-3., -3., 12.),
                RGB::new(0.8, 1., 1.),
                5.,
                0., 0., 0., 0.5,
                1.5
            )
        ]
    );

    let mut file = std::fs::File::create(format!("heckbert_{BOUNCES}_{VIEWING_ANGLE}.ppm")).unwrap();
    renderer.render_ppm(&mut file, &scene);
}