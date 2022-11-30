# Paul Heckbert's bussiness card renderer
Paul Heckbert's bussiness card ray tracer.  

## Package Files
The package consists of 3 files:

1. src/lib.rs  
Contains the implementation of Paul Heckbert's ray tracer. It's been heavily modified
and uses the external library `glm` for vectors.  
Used be the other executables.

2. examples/simple_scene.rs  
Builds the scene and calls the render method.  
Run using `cargo run --example simple_scene`, additionaly pass `--release` to run
with optimizations.

3. benches/depth-bench.rs  
Benchmarks the renderer with different depth values.  
Run using `cargo bench --bench depth-bench`.
