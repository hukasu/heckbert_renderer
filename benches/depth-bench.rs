use criterion::{Criterion, BenchmarkId};
use criterion::{criterion_main, criterion_group};
use heckbert_renderer::{Scene, Sphere, RGB, BussinessCardRenderer};

pub fn criterion_benchmark(c: &mut Criterion) {
    const SIZE: usize = 700;
    const VIEWING_ANGLE: f64 = 25.;

    let mut group = c.benchmark_group("Bussiness Card renderer");

    let spheres = vec![
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
    ];
    let scene = Scene::new(
        RGB::new(0.02, 0.02, 0.02),
        spheres
    );

    for depth in vec![1, 2, 3, 4, 5, 10, 20] {
        let id = BenchmarkId::new("Bussiness Card renderer", depth);
        group.bench_with_input(id, &depth, |b, depth| {
            b.iter(
                || {
                    let renderer = BussinessCardRenderer::new(
                        *depth,
                        (SIZE, SIZE),
                        VIEWING_ANGLE
                    );
                    renderer.render_ppm(&mut std::io::sink(), &scene);
                }
            )
        });
    }
    
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);