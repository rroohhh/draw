#[macro_use]
extern crate glium;
extern crate num;
extern crate rand;

use glium::{glutin, Surface};
use glutin::GlContext;


// TODO(robin): things to think about / decide
// - store brush in savefile 
// - use serde (serde bincode) for savefile?
// - store settings in savefile
// - how to handle bg color
// - calc bounding box for every stroke? (-> incremental bounding box)
// - level of detail for rendering speed improvement?
// - compression for savefile
//
// TODO(robin): features:
// - auto save
// - undo
// - undo tree?
// - layers?
// - markers?
// - smoothing

// 128 bits fixed point
// 4096 points -> 12 bits minimum
// -> 2^116 zoom

#[derive(Clone, Debug)]
struct CanvasPoint {
    x: i128,
    y: i128,
}

impl CanvasPoint {
    fn scale_by(&mut self, scale: &CanvasScale) -> &mut Self {
        self.x /= scale.x;      
        self.y /= scale.y;      

        self
    }
    
    fn offset_by(&mut self, offset: &CanvasPoint) -> &mut Self {
        self.x -= offset.x;
        self.y -= offset.y;

        self
    }
}

#[derive(Clone)]
struct CanvasScale {
    x: i128,
    y: i128,
}


impl CanvasPoint {
    fn to_nearest<T: num::Float>(&self) -> (T, T) {
        (num::NumCast::from(self.x).unwrap(), num::NumCast::from(self.y).unwrap())
    }
    
    fn to_nearest_with_offset_scale<T: num::Float>(&self, offset: &CanvasPoint, scale: &CanvasScale) -> (T, T) {
        // TODO(robin): think more about the order of this (they provide different guarantees
        // when limiting the minimum zoom to 12 bits this probably doesn't change much 
        let mut point = self.clone();
        point.offset_by(offset).scale_by(scale);

        (num::NumCast::from(point.x).unwrap(), num::NumCast::from(point.y).unwrap())
    }
}

#[derive(Clone)]
struct Color {}

#[derive(Clone)]
struct Tilt {
    x: f64, // normalized to [0.0, 2·π]
    y: f64, // normalized to [0.0, 2·π]
}

#[derive(Clone)]
struct StrokePoint {
    position: CanvasPoint,
    pressure: f64, // normalized to [0.0, 1.0]
    tilt: Tilt,
    rotation: f64, // normalized to [0.0, 2·π]
}

// one global array of CanvasPoints for all Strokes?
trait Stroke {
    fn add_points(&mut self, points: &[StrokePoint]);
    fn render(&self, context: &glium::backend::Facade, frame: &mut glium::Frame, offset: &CanvasPoint, scale: &CanvasScale);
    fn clear(&mut self);
    fn box_clone(&self) -> Box<Stroke>;
    fn empty(&self) -> bool;
}

#[derive(Clone)]
struct SimpleStroke {
    points: Vec<CanvasPoint>,
    pressure: Vec<f64>,
    color: Color,
    scale: CanvasScale
}

impl SimpleStroke {
    fn new(color: Color, scale: CanvasScale) -> SimpleStroke {
        SimpleStroke {
            points: Vec::new(),
            pressure: Vec::new(),
            color,
            scale,
        }
    }
}

impl Stroke for SimpleStroke {
    fn empty(&self) -> bool {
        self.points.len() == 0
    }

    fn box_clone(&self) -> Box<Stroke> {
        Box::new(self.clone())
    }

    fn clear(&mut self) {
        self.points.clear();
        self.pressure.clear();
    }

    fn add_points(&mut self, points: &[StrokePoint]) {
        for StrokePoint { position, pressure, .. } in points {
            self.points.push(position.clone());
            self.pressure.push(pressure.clone());
        }
    }

    fn render(&self, context: &glium::backend::Facade, frame: &mut glium::Frame, offset: &CanvasPoint, scale: &CanvasScale) {
        #[derive(Copy, Clone, Debug)]
        struct Vertex {
            position: [f32; 2],
            edge: f32,
        }

        // println!("{:?}", self.points);

        // println!("len: {}", self.points.len());
    
        implement_vertex!(Vertex, position, edge);
   
        // let line_points: Vec<_> = self.points.iter().map(|p| p.to_nearest_with_offset_scale::<f32>(offset, scale)).collect();
	let line_points: Vec<_> = self.points.iter().map(|p| (p.x, p.y)).collect();
        
        let mut points = Vec::new();

	let line_width = 4000.0;

 	if line_points.len() >= 2 {
            let (x1, y1) = line_points[0]; 
            let (x2, y2) = line_points[1]; 

            let dx = x2 - x1;
            let dy = y2 - y1;
            
            let len = (dx * dx + dy * dy);
            if len > 0 {
                let len = (len as f64).sqrt();

                let hat_dx = dx as f64 / len;
                let hat_dy = dy as f64 / len;

		let (x1, y1): (f64, f64) = self.points[0].to_nearest_with_offset_scale::<f64>(offset, scale);
                let dx = -hat_dy * line_width;
		let dy =  hat_dx * line_width;

                points.push(Vertex {
                    position: [(x1 - dx) as f32, (y1 - dy) as f32],
                    edge: -1.0
                });

                points.push(Vertex {
                    position: [(x1 - dx - dy) as f32, (y1 - dy + dx) as f32],
                    edge: -1.0
                });

                points.push(Vertex {
                    position: [x1 as f32, y1 as f32],
                    edge: 0.0
                });

                points.push(Vertex {
                    position: [(x1 + dx - dy) as f32, (y1 + dy + dx) as f32],
                    edge: -1.0
                });

                points.push(Vertex {
                    position: [(x1 + dx) as f32, (y1 + dy) as f32],
                    edge: -1.0
                });

                points.push(Vertex {
                    position: [(x1 + dx) as f32, (y1 + dy) as f32],
                    edge: -1.0
                });
                
                points.push(Vertex {
                    position: [(x1 - dx) as f32, (y1 - dy) as f32],
                    edge: 1.0
                });
            }
	}

        for i in 2..line_points.len() {
            let (x1, y1) = line_points[i - 2]; 
            let (x2, y2) = line_points[i - 1]; 
            let (x3, y3) = line_points[i]; 

            let dx_a = x2 - x1;
            let dy_a = y2 - y1;
            let len_a = ((dx_a * dx_a + dy_a * dy_a) as f64).sqrt();

            let dx_b = x3 - x2;
            let dy_b = y3 - y2;
            let len_b = ((dx_b * dx_b + dy_b * dy_b) as f64).sqrt();

            let dx = len_b * dx_a as f64 + len_a * dx_b as f64;
            let dy = len_b * dy_a as f64 + len_a * dy_b as f64;

            let len = (dx * dx + dy * dy);
            if len > 0.0 {
                let len = (len as f64).sqrt();
                let hat_dx = (dx as f64) / len;
                let hat_dy = (dy as f64) / len;

	        let factor = 1.0 / (hat_dx * dx_b as f64 + hat_dy * dy_b as f64) * len_b;
                let line_width = line_width * factor;
		
		if factor < 2.0 {
		    let (x1, y1): (f64, f64) = self.points[i - 2].to_nearest_with_offset_scale::<f64>(offset, scale);
		    let (x2, y2): (f64, f64) = self.points[i - 1].to_nearest_with_offset_scale::<f64>(offset, scale);
		    let (x3, y3): (f64, f64) = self.points[i].to_nearest_with_offset_scale::<f64>(offset, scale);
                    let xa = (x2 - hat_dy * line_width) as f64;
		    let ya = (y2 + hat_dx * line_width) as f64;

		    let xb = (x2 + hat_dy * line_width) as f64;
		    let yb = (y2 - hat_dx * line_width) as f64;

                    points.push(Vertex {
                        position: [xa as f32, ya as f32],
                        edge: -1.0
                    });
                    
                    points.push(Vertex {
                        position: [xb as f32, yb as f32],
                        edge: 1.0 
                    });
		} else {
		    // TODO(robin): make a flatted angle (or rounded?)
		    let (x1, y1): (f64, f64) = self.points[i - 2].to_nearest_with_offset_scale::<f64>(offset, scale);
		    let (x2, y2): (f64, f64) = self.points[i - 1].to_nearest_with_offset_scale::<f64>(offset, scale);
		    let (x3, y3): (f64, f64) = self.points[i].to_nearest_with_offset_scale::<f64>(offset, scale);

                    let xa = (x2 - hat_dy * line_width) as f64;
		    let ya = (y2 + hat_dx * line_width) as f64;

		    let xb = (x2 + hat_dy * line_width) as f64;
		    let yb = (y2 - hat_dx * line_width) as f64;

                    points.push(Vertex {
                        position: [xa as f32, ya as f32],
                        edge: -1.0
                    });
                    
                    points.push(Vertex {
                        position: [xb as f32, yb as f32],
                        edge: 1.0 
                    });

		}	
            }
        }

	if line_points.len() >= 2 {
            let i = line_points.len();
            let (x1, y1) = line_points[i - 2]; 
            let (x2, y2) = line_points[i - 1]; 

            let dx = x2 - x1;
            let dy = y2 - y1;

            let len = (dx * dx + dy * dy);
            if len > 0 {
                let len = (len as f64).sqrt();

                let hat_dx = (dx as f64) / len;
                let hat_dy = (dy as f64) / len;


		let (x2, y2): (f64, f64) = self.points[i - 1].to_nearest_with_offset_scale::<f64>(offset, scale);
                let dx = -hat_dy * line_width;
		let dy =  hat_dx * line_width;

                points.push(Vertex {
                    position: [(x2 + dx) as f32, (y2 + dy) as f32],
                    edge: -1.0
                });
                
                points.push(Vertex {
                    position: [(x2 - dx) as f32, (y2 + dy) as f32],
                    edge: 1.0
                });

                points.push(Vertex {
                    position: [(x2 - dx) as f32, (y2 - dy) as f32],
                    edge: -1.0
                });

                points.push(Vertex {
                    position: [(x2 - dx + dy) as f32, (y2 - dy - dx) as f32],
                    edge: -1.0
                });

                points.push(Vertex {
                    position: [x2 as f32, y2 as f32],
                    edge: 0.0
                });

                points.push(Vertex {
                    position: [(x2 + dx + dy) as f32, (y2 + dy - dx) as f32],
                    edge: -1.0
                });

                points.push(Vertex {
                    position: [(x2 + dx) as f32, (y2 + dy) as f32],
                    edge: -1.0
                });
            }
	}

        let vertex_buffer = glium::VertexBuffer::new(context, &points).unwrap();
        let indices = glium::index::NoIndices(glium::index::PrimitiveType::TriangleStrip);
    
        let vertex_shader_src = r#"
            #version 450
            in vec2 position;
            in float edge;
            out vec2 pos;
            out float edgedistance;

            void main() {
                edgedistance = edge;
                pos = position;
                gl_Position = vec4(position / 65536.0, 0.0, 1.0);
//                gl_PointSize = 4.0;
            }
        "#;
    
        let fragment_shader_src = r#"
            #version 450
            out vec4 color;
            in vec2 pos;
            in float edgedistance;

            void main() {
                vec2 pos = gl_PointCoord - 0.5;
                if (length(pos) > 0.8) {
                    discard;
                } else {
                    color = vec4(1.0, 0.0, 0.0, 1.0);
                }
                // color = vec4(1.0, 0.0, 0.0, mix(0.0, 1.0, 1.0 - abs(edgedistance)));
            }
        "#;
    	
        let program = glium::Program::new(context, glium::program::ProgramCreationInput::SourceCode {
            vertex_shader: vertex_shader_src,
            fragment_shader: fragment_shader_src,
            geometry_shader: None,
            tessellation_control_shader: None,
            tessellation_evaluation_shader: None,
            transform_feedback_varyings: None,
            uses_point_size: true, 
            outputs_srgb: false, 
        }).unwrap();
    
        let params = glium::DrawParameters {
            smooth: glium::draw_parameters::Smooth::Nicest.into(),
            line_width: 5.0.into(),
            blend: glium::draw_parameters::Blend::alpha_blending(),
            .. Default::default()
        };
        
        frame.draw(&vertex_buffer, &indices, &program, &glium::uniforms::EmptyUniforms,
                    &params).unwrap();
    }
}

fn main() {
    let mut events_loop = glutin::EventsLoop::new();
    let window = glutin::WindowBuilder::new().with_title("draw");
    let context = glutin::ContextBuilder::new().with_multisampling(4);
    let context = glutin::ContextBuilder::new();
    let display = glium::Display::new(window, context, &events_loop).unwrap();
    
    let mut dpi_factor = display.gl_window().window().get_hidpi_factor();
    let mut closed = false;

    let mut stroke = SimpleStroke::new(Color {}, CanvasScale { x: 1, y: 1 });

    use rand::{thread_rng, Rng};

    /*
    for i in 0..1000 {
        let x: i64 = thread_rng().gen_range(-8, 8);
        let y: i64 = thread_rng().gen_range(-8, 8);

        stroke.add_points(&[ StrokePoint { 
            position: CanvasPoint { x: x as i128, y: y as i128 },
            pressure: 0.0,
            tilt: Tilt { x: 0.0, y: 0.0 },
            rotation: 0.0,
        } ]);
    }
    */

    let mut mouse_down = false;
    let mut strokes: Vec<Box<Stroke>> = Vec::new();


// TODO(robin): figure out the overlay situation 
// test case:
// [CanvasPoint { x: 10, y: -7 }, CanvasPoint { x: 9, y: -6 }, CanvasPoint { x: 7, y: -4 }, CanvasPoint { x: 4, y: -2 }, CanvasPoint { x: 1, y: 1 }, CanvasPoint { x: 0, y: 2 }, CanvasPoint { x: -5, y: 5 }, CanvasPoint { x: -5, y: 6 }, CanvasPoint { x: -5, y: 5 }, CanvasPoint { x: -6, y: 5 }, CanvasPoint { x: -6, y: 3 }, CanvasPoint { x: -7, y: 2 }, CanvasPoint { x: -7, y: 1 }, CanvasPoint { x: -7, y: 0 }]

/*
        let canvas_points =  [CanvasPoint { x: -1, y: -4 }, CanvasPoint { x: -2, y: -6 }, CanvasPoint { x: -2, y: -7 }, CanvasPoint { x: -3, y: -8 }, CanvasPoint { x: -4, y: -8 }, CanvasPoint { x: -4, y: -7 }];
        for point in &canvas_points[3..] {
            stroke.add_points(&[ StrokePoint {
                position: point.clone(),
                pressure: 0.0, 
                tilt: Tilt { x: 0.0, y: 0.0 },
                rotation: 0.0,
            } ]);
	}
*/



/*
                            stroke.add_points(&[ StrokePoint {
                                position: CanvasPoint { x: -5 as i128, y: 0 as i128 },
                                pressure: 0.0, 
                                tilt: Tilt { x: 0.0, y: 0.0 },
                                rotation: 0.0,
                            } ]);

                            stroke.add_points(&[ StrokePoint {
                                position: CanvasPoint { x: 0 as i128, y: 0 as i128 },
                                pressure: 0.0, 
                                tilt: Tilt { x: 0.0, y: 0.0 },
                                rotation: 0.0,
                            } ]);

                            stroke.add_points(&[ StrokePoint {
                                position: CanvasPoint { x: 0 as i128, y: 5 as i128 },
                                pressure: 0.0, 
                                tilt: Tilt { x: 0.0, y: 0.0 },
                                rotation: 0.0,
                            } ]);
*/

//    strokes.push(stroke.box_clone());
//
    let mut last_pos = (-10923874109324719820374, -2134812094312349);

    while !closed {
        let mut target = display.draw();
        target.clear_color(0.0, 0.0, 1.0, 1.0);
        for stroke in strokes.iter() {
            stroke.render(&display, &mut target, &CanvasPoint { x: 0, y: 0 }, &CanvasScale { x: 1, y: 1 });
        }

        stroke.render(&display, &mut target, &CanvasPoint { x: 0, y: 0 }, &CanvasScale { x: 1, y: 1 });
        target.finish().unwrap();


        events_loop.poll_events(|event| {
            match event {
                glutin::Event::WindowEvent { event, .. } => match event {
                    glutin::WindowEvent::HiDpiFactorChanged(factor) => {
                        dpi_factor = factor;
                    }
                    glutin::WindowEvent::CursorMoved { device_id, position, modifiers } => {
                        // println!("cursor moved {:?} {:?} {:?}", device_id, position, modifiers);
                         // self.display.get_framebuffer_dimensions
                         let (dimx, dimy) = display.get_framebuffer_dimensions();
                         let dimx = dimx as f64;
                         let dimy = dimy as f64;
                         let position = position.to_physical(dpi_factor);
                         // println!("mouse {} {}", position.x, position.y);
                         // println!("dim   {} {}", dimx, dimy);
                         // println!("dimx {}, dimy {}", dimx, dimy);
                         let mut x = (position.x / dimx) * 2.0 - 1.0;
                         let mut y = ((dimy - position.y) / dimy) * 2.0 - 1.0;
                         x *= 65536.0;
                         y *= 65536.0;
                         let x = x.round() as i128;
                         let y = y.round() as i128;


                         if mouse_down {
                            if (x, y) != last_pos {
                                last_pos = (x, y);
                                stroke.add_points(&[ StrokePoint {
                                    position: CanvasPoint { x: x as i128, y: y as i128 },
                                    pressure: 0.0, 
                                    tilt: Tilt { x: 0.0, y: 0.0 },
                                    rotation: 0.0,
                                } ]);
                            }
                         }
                    }
                    glutin::WindowEvent::MouseInput { device_id, state, button, modifiers } => {
                        if button == glium::glutin::MouseButton::Left {
                            if state == glium::glutin::ElementState::Pressed {
                                mouse_down = true;
                            } else {
                                mouse_down = false;
                                if !stroke.empty() {
                                    strokes.push(stroke.box_clone());
                                    stroke.clear();
                                }
                            }
                        }
                    }
                    glutin::WindowEvent::AxisMotion { device_id, axis, value } => {
                        if axis == 2 {
                            if value > 1000.0 {
                                mouse_down = true;
                            } else {
                                mouse_down = false;
                                if !stroke.empty() {
                                    strokes.push(stroke.box_clone());
                                    stroke.clear();
                                }
                            }
                        }
//                          println!("axis motion {:?} {:?} {:?}", device_id, axis, value);
                    }
                    glutin::WindowEvent::CloseRequested => closed = true,
                    _ => ()
                },
                _ => (),
            }
        });
    }
}
