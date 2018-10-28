# draw
nearly infinite canvas vector drawing program, mainly for keeping a digital notebook

## technology
every point drawn is represented by a 128bit integer, this gives a dynamic range of over 3*10^38.
using integer means that accuracy is not reduced when moving away from the origin (you might recognise this behaviour from other "infinite" canvas drawing programs)

points belonging to a stroke are connected with triangle strips and the stroke shape is rendered in a pixel shader

## goals
- high performance: using this for normal notetaking no slowdowns should be noticed
- collaborative: many people should be able to draw together transparently over network
- multiplatform: native linux, windows, mac, android and webbrowser support
- infinite undo / save: don't worry about saving 
- full featured: all the normal drawing tools are available: curves, straight lines, polygons, layers, markers, rotation, translation
- good tablet support: support drawing tablets with pressure, tilt and rotation
- flexible output: high quality svg and png rendering 
