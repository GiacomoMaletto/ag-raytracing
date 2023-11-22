local shader = love.graphics.newShader([[
    uniform vec3 pos;
    uniform vec3 dir;
    uniform vec3 x_dir;
    uniform vec3 y_dir;

    float ev(vec4 a, float t) {
        return ((a[3]*t + a[2])*t + a[1])*t + a[0];
    }

    int var(vec4 a) {
        vec4 as = vec4(a[3] + a[2] + a[1] + a[0], a[2] + 2.0*a[1] + 3.0*a[0], a[1] + 3.0*a[0], a[0]);
        int s = 0;
        int v = 0;
        for (int i = 0; i < 4; i++) {
            if (as[i] > 0.0) {
                if (s < 0) v++;
                s = 1;
            } else if (as[i] < 0.0) {
                if (s > 0) v++;
                s = -1;
            }
        }
        return v;
    }

    vec4 effect(vec4 color, Image texture, vec2 tc, vec2 screen_coords){
        float t = 0.0;

        float x = 2.0*tc.x - 1.0;
        float y = 1.0 - 2.0*tc.y;
        vec3 v = normalize(dir + x*x_dir + y*y_dir);
        // vec3 v = dir + x*x_dir + y*y_dir;

        float x0 = pos.x;
        float y0 = pos.y;
        float z0 = pos.z;

        float X = v.x;
        float Y = v.y;
        float Z = v.z;

        // vec4 a = vec4(x0*y0*z0 + x0 + 1.0,
        //               Z*x0*y0 + Y*x0*z0 + X*y0*z0 + X,
        //               Y*Z*x0 + X*Z*y0 + X*Y*z0,
        //               X*Y*Z);

        int degree = 3;

        // vec4 a = vec4(x0*x0 + y0*y0 + z0*z0 - 1,
        //               2*(X*x0 + Y*y0 + Z*z0),
        //               X*X + Y*Y + Z*Z,
        //               0);

        vec4 a = vec4(+ (-2*x0*y0*y0 + 3*x0*x0*z0 + 2*y0*y0*z0 - 2*x0*z0*z0 + x0*y0 - y0*z0 - z0),
                      + (3*Z*x0*x0 - 4*Y*x0*y0 - 2*X*y0*y0 + 2*Z*y0*y0 + 6*X*x0*z0 - 4*Z*x0*z0 + 4*Y*y0*z0 - 2*X*z0*z0 + Y*x0 + X*y0 - Z*y0 - Y*z0 - Z),
                      - (2*Y*Y*x0 - 6*X*Z*x0 + 2*Z*Z*x0 + 4*X*Y*y0 - 4*Y*Z*y0 - 3*X*X*z0 - 2*Y*Y*z0 + 4*X*Z*z0 - X*Y + Y*Z),
                      -(2*X*Y*Y - 3*X*X*Z - 2*Y*Y*Z + 2*X*Z*Z));

        vec4 ap = vec4(a[1], 2.0*a[2], 3.0*a[3], 0.0);

        float b = 0.0;
        for (int i = 0; i < degree; i++) {
            float ratio = -a[i]/a[degree];
            if (ratio > 0) b = max(b, pow(ratio, 1.0/(3 - i)));
        }
        b *= 2.0;

        vec4 as = vec4(a[0], b*a[1], b*b*a[2], b*b*b*a[3]);

        vec4 P[20];
        float S[20];
        float E[20];
        P[0] = as;
        S[0] = 0.0;
        E[0] = 1.0;
        int last = 0;
        int n = 0;

        vec4 A;
        float s, e;

        while (n < 20 && last >= 0 && last < 20 - 2) {
            A = P[last];
            s = S[last];
            e = E[last];

            int V = var(A);
            if (V == 0) {
                last--;
                if (abs(ev(as, e)) < 0.00001) { t = e*b; break; }
            } else {
                t = s*b;

                vec4 Ap = vec4(A[0], A[1] / 2.0, A[2] / 4.0, A[3] / 8.0);
                // vec4 Ap = vec4(A[0] * 8.0, A[1] * 4.0, A[2] * 2.0, A[3]);
                P[last + 1] = Ap;
                S[last + 1] = s;
                E[last + 1] = (s + e)/2.0;

                P[last] = vec4(Ap[0] + Ap[1] + Ap[2] + Ap[3], Ap[1] + 2.0*Ap[2] + 3.0*Ap[3], Ap[2] + 3.0*Ap[3], Ap[3]);
                S[last] = (s + e)/2.0;
                E[last] = e;

                last++;
            }

            n++;
        }

        if (last < 0) t = 0.0;
        else {
            for (int i = 0; i < 10; i++) {
                t -= ev(a, t)/ev(ap, t);
                if (t < s*b || t > e*b) { t = s*b; break; }
            }
        }

        vec4 type;
        if (n >= 30) type = vec4(1.0, 0.0, 0.0, 1.0);
        else if (last < 0) type = vec4(0.0, 1.0, 0.0, 1.0);
        else if (last >= 30 - 2) type = vec4(0.0, 0.0, 1.0, 1.0);
        else type = vec4(1.0, 1.0, 1.0, 1.0);

        float c;
        if (t > 0) {
            // c = pow(max(0.1*t, 0.0) + 1.0, -2.0);
            c = exp(-0.2*t);
            // c = pow(max(0.3*t, 0.0) + 1.0, -1.0);
        } else { c = 0.0; }

        return (vec4(c, c, c, 1.0));
    }
]])

local V, M = unpack(require "vector")

love.mouse.setRelativeMode(true)
love.graphics.setDefaultFilter("nearest", "nearest")

local sw, sh = love.graphics.getDimensions()

local camera = {}
camera.xyz = {5, 5, 5}
camera.phi = math.rad(180 + 45)
camera.theta = math.rad(120)
camera.speed = 1
camera.fov = math.rad(120)
function camera.dir()
    return V.spherical(1, camera.theta, camera.phi)
end
function camera.right()
    return V.scale(math.tan(camera.fov/2), V.unit(V.cross(camera.dir(), {0,0,1})))
end
function camera.up()
    return V.scale(sh/sw, V.cross(camera.right(), camera.dir()))
end

local white_canvas = love.graphics.newCanvas(sw, sh, {format="rgba32f"})
love.graphics.setCanvas(white)
love.graphics.clear(1, 1, 1)
love.graphics.setCanvas()

local distance_canvas = love.graphics.newCanvas(sw, sh, {format="rgba32f"})

function love.mousemoved(x, y, dx, dy, istouch)
    camera.phi = camera.phi - dx/200
    local new_theta = camera.theta + dy/200
    if 0 <= new_theta and new_theta <= math.rad(180) then camera.theta = new_theta end
end

local dt = 1/60
local t = 0
function love.update(Dt)
    dt = Dt
    t = t + dt
    if love.keyboard.isDown("escape") then
        love.event.quit()
    end

    if love.keyboard.isDown("w") then
        camera.xyz = V.add(camera.xyz, V.mul(dt*camera.speed, camera.dir()))
    end
    if love.keyboard.isDown("s") then
        camera.xyz = V.add(camera.xyz, V.mul(-dt*camera.speed, camera.dir()))
    end
    if love.keyboard.isDown("d") then
        camera.xyz = V.add(camera.xyz, V.mul(dt*camera.speed, camera.right()))
    end
    if love.keyboard.isDown("a") then
        camera.xyz = V.add(camera.xyz, V.mul(-dt*camera.speed, camera.right()))
    end
    if love.keyboard.isDown("e") then
        camera.xyz = V.add(camera.xyz, V.mul(dt*camera.speed, camera.up()))
    end
    if love.keyboard.isDown("q") then
        camera.xyz = V.add(camera.xyz, V.mul(-dt*camera.speed, camera.up()))
    end

    shader:send("pos", camera.xyz)
    shader:send("dir", camera.dir())
    shader:send("x_dir", camera.right())
    shader:send("y_dir", camera.up())
end

function love.draw()
    love.graphics.setShader(shader)
    love.graphics.setCanvas(distance_canvas)
    love.graphics.draw(white_canvas)
    love.graphics.setCanvas()
    love.graphics.setShader()

    love.graphics.draw(distance_canvas)
end