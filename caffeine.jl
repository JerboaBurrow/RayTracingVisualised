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
Draw a sphere
"""
function drawSphere!(ax, centre::Point3f, radius::Float32, colour::Observable{RGBAf})
    U = LinRange(0.0, pi, 100)
    V = LinRange(0.0, 2*pi, 100)

    x = [radius*sin(u)*cos(v) for u in U, v in V].+centre[1]
    y = [radius*sin(u)*sin(v) for u in U, v in V].+centre[2]
    z = [radius*cos(u) for u in U, v in V].+centre[3]
    surface!(ax, x, y, z; color=colour, shading = NoShading, transparency = true, alpha=0.125)
    wireframe!(ax, x, y, z; color=colour, linewidth = 0.5, transparency = true, alpha=0.125)
end

"""
Draw a pixel grid out of oriented quads. Pixels are coloured based
  upon simple ray-tracing of a collection of spheres.
"""
function pixelGrid!(ax, spherePositions::Vector{Point3f}, sphereRadii::Vector{Float32}, sphereColours::Vector{Observable{RGBAf}}, rayOrigin::Point3f; gridOffset::Float32=Float32(1.0), gridPoints::Int64=16, gridSize::Float32=Float32(0.1))
    rayDirection = -rayOrigin
    rayDirection/=norm(rayDirection)
    gridOrigin = rayOrigin+rayDirection*gridOffset
    base, q_a, q_b = quadVertices(rayDirection, gridOrigin, Float32(pi/4.0), gridSize)
    q_a = Point3f(q_a*gridSize)
    q_b = Point3f(q_b*gridSize)
    colours = [Observable(RGBAf(0.5,0.5,0.5,0.0)) for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    rayColours = [Observable(RGBAf(0.0,0.0,0.0,0.0)) for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    rayDirections = [Point3f(0) for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    rayLighting = [0.0 for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    rayIntersections = [((Point3f(NaN), Point3f(NaN)), (Point3f(NaN), Point3f(NaN))) for i in -gridPoints+1:gridPoints-1 for j in -gridPoints+1:gridPoints-1]
    k=1
    for i in -gridPoints+1:gridPoints-1
        for j in -gridPoints+1:gridPoints-1
            pos = gridOrigin+q_a*i+q_b*j
            r = pos-rayOrigin
            r/=norm(r)
            rayDirections[k] = r
            for (spherePos, radius, sphereColour) in zip(spherePositions, sphereRadii, sphereColours)
                (ai, an), (bi, bn) = intersections(spherePos, radius, r, rayOrigin)
                if !isnan(ai)
                    l = rayOrigin-spherePos
                    l /= sqrt(sum(l.*l))
                    rayLighting[k] = 0.0001+0.75*max(sum(an.*l), 0.0)
                    colours[k] = Observable(RGBAf(sphereColour[].r, sphereColour[].g, sphereColour[].b, 0.0))
                    rayIntersections[k] = ((ai, an), (bi, bn))
                    break
                end
            end
            quad!(ax, rayDirection, position=pos, theta=Float32(pi/4.0), scale=gridSize, colour=colours[k])
            lines!(ax, [rayOrigin, rayOrigin+r*4], color=rayColours[k])
            k+=1
        end
    end
    return colours, rayColours, rayDirections, rayLighting, rayIntersections
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

# Camera position.
cz = 6.0
p0 = [cz, cz, cz]

# The molecule's data.
caf = []
radii = Dict([("H", Float32(1.2)), ("N", Float32(1.66)), ("C", Float32(1.77)), ("O", Float32(1.5))])
colours = Dict([("H", RGBAf(1.0,1.0,1.0,1.0)), ("N", RGBAf(0.561,0.561,1.0,1.0)), ("C", RGBAf(0.784,0.784,0.784,1.0)), ("O", RGBAf(0.941,0.0,0.0,1.0))])
for line in readlines("caffiene.xyz")[3:end]
    r = radii[split(line)[1]]
    colour = colours[split(line)[1]]
    centre = Point3f(parse.(Float32, (String.(split(line)[2:end]))))
    push!(caf, [centre, r, Observable(colour)])
end

# Centre it.
com = Point3f(0)
for c in caf com += c[1] end
com /= length(caf)
for i in 1:length(caf) caf[i][1]-=com end

# Distances for render order.
dists = [norm(c[1]-p0) for c in caf];
indices = sortperm(dists)
caf = caf[indices]

positions = [c[1] for c in caf]
radii = [c[2] for c in caf]
colours = [c[3] for c in caf]

[drawSphere!(ax, c[1], c[2]*Float32(0.5), c[3]) for c in caf]

pixelColours, rayColours, rays, rayLighting, rayIntersections = pixelGrid!(ax, positions, radii, colours, Point3f(p0), gridPoints=25, gridSize=Float32(0.025))

scatter!(ax, Point3f(p0), markersize=10, color = :red)

a = Observable(Point3f(NaN))
b = Observable(Point3f(NaN))
scatter!(ax, a, markersize=10, color = :red)
scatter!(ax, b, markersize=10, color=:red)

p0 = Point3f(p0)
p1 = Point3f(p0)
ps = Observable([p0, p1])
lines!(ax, ps, color=:red)

normal = Observable([p0, p1])
lines!(ax, normal, color=:white)

L = cz+1
xlims!(ax, -L, L)
ylims!(ax, -L, L)
zlims!(ax, -L, L)

fps=60

queue = collect(1:length(rays))

d = norm(p0)*2.0
currentRay = pop!(queue)
n = rays[currentRay]

cam = cam3d!(ax.scene)
rayTime = 0.1
t = 0.0
intersection = false

hidedecorations!(ax)

times = range(0.0, rayTime*(length(rays)+50), step=1.0/fps)
dt = 1.0/fps
record(fig, "ray.mp4", times; framerate = fps) do time
    if time == 0.0
        rotate_cam!(ax.scene, Vec3f(0.0, 2.0*pi*1.0, 0.0))
    end
    ps[] = [p0, ps[][1]+Float32(t*d/rayTime)*n]
    global t += dt

    if t > rayTime*0.5 && !intersection
        (ai, an), (bi, bn) = rayIntersections[currentRay]
        a[] = ai
        b[] = bi
        normal[] = [ai, an*1.25]
        global intersection = true
    end

    if t > rayTime && length(queue) > 0
        print(length(queue)," ")
        colour = pixelColours[currentRay][]
        if (rayLighting[currentRay] > 0.0)
            pixelColours[currentRay][] = RGBAf(rayLighting[currentRay]*colour.r, rayLighting[currentRay]*colour.g, rayLighting[currentRay]*colour.b, 1.0)
        else
            pixelColours[currentRay][] = RGBAf(colour.r, colour.g, colour.b, 1.0)
        end
        a[] = Point3f(NaN)
        b[] = Point3f(NaN)
        normal[] = [Point3f(NaN), Point3f(NaN)]
        global currentRay = pop!(queue)
        global n = rays[currentRay]
        global t = 0.0
        global intersection = false
    end
end