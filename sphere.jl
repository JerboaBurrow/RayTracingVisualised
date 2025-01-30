# Renders a visualisation of a sphere being ray-traced.
# A pixel grid is filled in by cast rays over time.

using GLMakie, GeometryBasics, LinearAlgebra, Random
# CPU will take hours rather than minutes.
# using CairoMakie

"""
Determine the intersection points on a sphere of a ray from the
  camera.

Returns also the normals vectors.
"""
function intersections(c, r, ray, camera)
    b = sum(camera.*ray)-sum(ray.*c)
    r2 = r*r
    det = b*b-(sum(c.*c)+sum(camera.*camera)-r2-2.0*sum(c.*camera))
    if det < 0
        return (Point3f(NaN,NaN,NaN), Point3f(NaN,NaN,NaN)), (Point3f(NaN,NaN,NaN), Point3f(NaN,NaN,NaN))
    end
    det = sqrt(det)
    a = ray * (-b+det)
    b = ray * (-b-det)
    an = a-c+camera
    bn = b-c+camera
    an = an/sqrt(sum(an.*an))
    bn = bn/sqrt(sum(bn.*bn))
    return (a+camera, an), (b+camera, bn)
end

"""
Determine the rotation matrix rotating a vector, from, to align with
  an axis, to.
"""
function rotateTo(from, to)
    # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    v = cross(from, to)
    vx = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
    return I(3) + vx + (vx*vx)/(1.0+sum(from.*to))
end

"""
Vertices of an oriented quad with normal vector normal and
  rotation theta about the normal axis.
"""
function quadVertices(normal::Point3f, position::Point3f, theta::Float32 = Float32(0.0), scale::Float32 = Float32(1.0))
    # Quad faces normal
    c = cos(theta)
    s = sin(theta)
    R = rotateTo(Point3f(0, 0, 1), normal)*[c -s 0; s c 0; 0 0 1]
    a = R*Point3f(1,0,0)
    b = R*Point3f(0,1,0)
    # Rotate, scale, translate
    return [Point3f(R*point)*scale+position for point in [Point3f(-0.5, -0.5, 0), Point3f(0.5, -0.5, 0), Point3f(0.5, 0.5, 0), Point3f(-0.5, 0.5, 0)]], a, b
end

"""
Draw an oriented quad.
"""
function quad!(ax, normal::Point3f; position::Point3f, theta::Float32 = Float32(0.0), scale::Float32 = Float32(1.0), colour::Observable{RGBAf} = RGBAf(1.0, 0.75, 0.75,1.0))
    q, = quadVertices(normal, position, theta, scale)
    uv_mesh = map((p) -> GeometryBasics.Mesh(p, GLTriangleFace[(1, 2, 3), (3, 4, 1)]; uv=Vec3f.(p)), Observable(q))
    mesh!(ax, uv_mesh, color=colour, shading=NoShading)
end

"""
Draw a pixel grid out of oriented quads. Pixels are coloured based
  upon simple ray-tracing of a sphere.
"""
function pixelGrid!(ax, spherePos::Point3f, rayOrigin::Point3f; radius::Float32=Float32(1.0), gridOffset::Float32=Float32(1.0), gridPoints::Int64=16, gridSize::Float32=Float32(0.1))
    rayDirection = spherePos-rayOrigin
    rayDirection/=norm(rayDirection)
    gridOrigin = rayOrigin+rayDirection*gridOffset
    base, q_a, q_b = quadVertices(rayDirection, gridOrigin, Float32(pi/4.0), gridSize)
    q_a = Point3f(q_a*gridSize)
    q_b = Point3f(q_b*gridSize)
    colours = [Observable(RGBAf(1.0,1.0,1.0,0.0)) for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    rayColours = [Observable(RGBAf(0.0,0.0,0.0,0.0)) for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    rayDirections = [Point3f(0) for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    k=1
    for i in -gridPoints+1:gridPoints-1
        for j in -gridPoints+1:gridPoints-1
            pos = gridOrigin+q_a*i+q_b*j
            r = pos-rayOrigin
            r/=norm(r)
            rayDirections[k] = r
            (ai, an), (bi, bn) = intersections(spherePos, radius, r, rayOrigin)
            if !isnan(ai)
                colours[k][] = RGBAf(1.0, 0.0, 0.0, 0.0)
            end
            quad!(ax, rayDirection, position=pos, theta=Float32(pi/4.0), scale=gridSize, colour=colours[k])
            lines!(ax, [rayOrigin, rayOrigin+r*4], color=rayColours[k])
            k+=1
        end
    end
    return colours, rayColours, rayDirections
end

set_theme!(backgroundcolor = :black)
fig = Figure(size=(1920,1080),fxaa=true, figure_padding=1)
ax = Axis3(
    fig[1, 1],
    aspect = :data,
    perspectiveness = 0.5,
    elevation = π / 9,
    yspinecolor_1=:white,
    yspinecolor_2=:white,
    yspinecolor_3=:transparent,
    yspinecolor_4=:transparent,
    xspinecolor_4=:transparent,
    xspinecolor_3=:transparent,
    xspinecolor_2=:white,
    xspinecolor_1=:white,
    zspinecolor_1=:transparent,
    zspinecolor_2=:white,
    zspinecolor_3=:transparent,
    zspinecolor_4=:white,
    xspinewidth=2,
    yspinewidth=2,
    zspinewidth=2
)

# A sphere mesh.
r = 1
U = LinRange(0.0, pi, 100)
V = LinRange(0.0, 2*pi, 100)

x = [r*sin(u)*cos(v) for u in U, v in V]
y = [r*sin(u)*sin(v) for u in U, v in V]
z = [r*cos(u) for u in U, v in V];

# Camera position.
cz = 2.0
p0 = [cz, cz, cz]

# Visualisation of a ray tracers result.
pixelColours, rayColours, rays = pixelGrid!(ax, Point3f(0), Point3f(p0), gridPoints=8, gridSize=Float32(0.05), radius=Float32(r))

# Camera point.
scatter!(ax, Point3f(p0), markersize=10, color = :red)

# Intersection points.
a = Observable(Point3f(NaN))
b = Observable(Point3f(NaN))
scatter!(ax, a, markersize=10, color = :red)
scatter!(ax, b, markersize=10, color=:red)

# Rays.
p0 = Point3f(p0)
p1 = Point3f(p0)
ps = Observable([p0, p1])
lines!(ax, ps, color=:red)

# To keep it drawn last.
x = Observable(x)
# Sphere and frame.
surface!(ax, x, y, z; colormap = :viridis, shading = NoShading,
    transparency = true, alpha=0.125)
wireframe!(ax, x, y, z; linewidth = 0.5, transparency = true, alpha=0.125)

xlims!(ax, -cz, cz)
ylims!(ax, -cz, cz)
zlims!(ax, -cz, cz)

fig[1,1] = ax

fps=60

# Each ray is drawn one at a time.
queue = randperm(length(rays))

d = norm(p0)*2.0
currentRay = pop!(queue)
n = rays[currentRay]

cam = cam3d!(ax.scene)
rayTime = 2.0
t = 0.0
intersection = false

hidedecorations!(ax)

times = range(0.0, rayTime*(length(rays)+5), step=1.0/fps)
dt = times[2]-times[1]
record(fig, "ray.mp4", times; framerate = fps) do time
    if time == 0.0
        rotate_cam!(ax.scene, Vec3f(0.0, 2.0*pi*1.05, 0.0))
    end
    rotate_cam!(ax.scene, Vec3f(0.0, dt*4.0*2.0*pi/60.0, 0.0))
    ps[] = [p0, ps[][1]+Float32(t*d/rayTime)*n]
    global t += dt

    if t > rayTime*0.5 && !intersection
        (ai, an), (bi, bn) = intersections(Point3f(0), Float32(r), rays[currentRay], p0)
        a[] = ai
        b[] = bi
        normal[] = [ai, an*1.25]
        global intersection = true
    end

    if t > rayTime && length(queue) > 0
        print(length(queue))
        colour = pixelColours[currentRay][]
        pixelColours[currentRay][] = RGBAf(colour.r, colour.g, colour.b, 1.0)
        a[] = Point3f(NaN)
        b[] = Point3f(NaN)
        normal[] = [Point3f(NaN), Point3f(NaN)]
        global currentRay = pop!(queue)
        global n = rays[currentRay]
        global t = 0.0
        global intersection = false
    end
    x[] = x[]
end
