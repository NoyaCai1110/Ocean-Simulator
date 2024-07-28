#version 330 core
 
in vec3 ourColor;
in vec2 screenCoord;
out vec4 FragColor;
uniform float timeStamp;
uniform float time_t;
uniform vec3 mousePos;
float DayTime = time_t/10;

float PAI = 3.14159265357989;
vec3 WaterCenter = vec3(0.0, -0.4, -2.0);

struct Light
{
    vec3 position;
    vec3 color;
};

struct Ray 
{
    vec3 origin;
    vec3 direction;
}; 

struct Sphere 
{
    vec3 center;
    float radius;
    vec3 amb_c;
    vec3 dif_c;
    vec3 spec_c;
    float ka;   //Ambient
    float kd;   //Diffuse
    float ks;   //Specular
    int ke;   //Specular coef
    float kr;   //Reflection
    float kt;   //Transmission
    float eta;  //Transmission refractive index
}; 

struct Triangle
{
    vec3 v1;
    vec3 v2;
    vec3 v3;
    vec3 amb_c;
    vec3 dif_c;
    vec3 spec_c;
    float ka;   //Ambient
    float kd;   //Diffuse
    float ks;   //Specular
    int ke;     //Specular coef
    float kr;   //Reflection
    float kt;   //Transmission
    float eta;  //Transmission refractive index
    bool iswave;//surface is wave
    float kw;   //determined by wave postion and how it affect the luminance
};

Light initLight(vec3 pos, vec3 color)
{
    Light light;
    light.position = pos;
    light.color = color;
    return light;
}

Ray initRay(vec3 origin, vec3 direction)
{
	Ray ray;
	ray.origin = origin;
	ray.direction = normalize(direction);

	return ray;
}

vec3 ExtendedRay(Ray ray, float t)
{
	return ray.origin + t * ray.direction;
}

Sphere initSphere(vec3 center, float radius, vec3 amb_c, vec3 dif_c, vec3 spec_c, float ka, float kd, float ks, int ke, float kr, float kt, float eta)
{
    Sphere sphere;
    sphere.center = center;
    sphere.radius = radius;
    sphere.amb_c = amb_c;
    sphere.dif_c = dif_c;
    sphere.spec_c = spec_c;
    sphere.ka = ka;   //Ambient
    sphere.kd = kd;   //Diffuse
    sphere.ks = ks;   //Specular
    sphere.ke = ke;   //Specular coef
    sphere.kr = kr;   //Reflection   
    sphere.kt = kt;   //Transmission
    sphere.eta = eta; //Transmission refractive index
    return sphere;
}

Triangle initTriangle(vec3 v1, vec3 v2, vec3 v3, vec3 amb_c, vec3 dif_c, vec3 spec_c, float ka, float kd, float ks, int ke, float kr, float kt, float eta, bool iswave)
{
    Triangle triangle;
    triangle.v1 = v1;
    triangle.v2 = v2;
    triangle.v3 = v3;
    triangle.amb_c = amb_c;
    triangle.dif_c = dif_c;
    triangle.spec_c = spec_c;
    triangle.ka = ka;   //Ambient
    triangle.kd = kd;   //Diffuse
    triangle.ks = ks;   //Specular
    triangle.ke = ke;   //Specular coef
    triangle.kr = kr;   //Reflection   
    triangle.kt = kt;   //Transmission
    triangle.eta = eta; //Transmission refractive index
    triangle.iswave = iswave;
    triangle.kw = 1;
    return triangle;
}

vec3 GerstnerWave (vec4 wave, vec3 p, inout vec3 tangent, inout vec3 binormal) 
{
    float steepness = wave.z;
    float wavelength = wave.w;
    float k = 2 * PAI / wavelength;
    float c = sqrt(9.8 / k);                // _WaveSpeed
    vec2 dir = normalize(wave.xy);
    float f = k * (dot(dir, p.xz) - c * time_t/3);
    float a = steepness / k;                // _Amplitude

    tangent += vec3(-dir.x * dir.x * (steepness * sin(f)), dir.x * (steepness * cos(f)), -dir.x * dir.y * (steepness * sin(f)));
    binormal += vec3(-dir.x * dir.y * (steepness * sin(f)), dir.y * (steepness * cos(f)), -dir.y * dir.y * (steepness * sin(f)));

    return vec3(dir.x * (a * cos(f)), a * sin(f), dir.y * (a * cos(f)));
}

/*Light sources and Primitive Shapes Init BEGIN*/
vec3 MYcolor = vec3((sin(time_t + 0 * PAI / 3) + 1) / 2, (sin(time_t + 2 * PAI / 3) + 1) / 2, (sin(time_t + 4 * PAI / 3) + 1) / 2);
float sin_t = (sin(time_t/3) + 1) / 2;
vec3 Skycolor = vec3(0.6, 0.8, 0.95);
vec3 SunPos = vec3(18*cos(DayTime), 18.0*sin(DayTime), -20.0);
vec3 MoonPos = vec3(-18*cos(DayTime), -18.0*sin(DayTime), -20.0);
float SunH = 18.0*sin(DayTime);
vec3 CurPos = vec3(18*cos(DayTime)* (SunH+0.4)/abs(SunH+0.4), 18.0*sin(DayTime) * (SunH+0.4)/abs(SunH+0.4), -20.0);
Light light1 = initLight(vec3(0.0, 5.0, 1.0), vec3(1.0, 1.0, 1.0));
Light light2 = initLight(vec3(2.0, 5.0, -1.0), vec3(1.0, 1.0, 1.0));
Light Sun = initLight(SunPos, vec3(1.0, 1.0, 1.0));
Light Moon = initLight(MoonPos, vec3(1.0, 1.0, 1.0));
Light Cur_light = initLight(CurPos, vec3(1.0, 1.0, 1.0));
//Light Cur_light = initLight(vec3(0.0, 18.0, -20.0), vec3(1.0, 1.0, 1.0));
Light l_list[1] = Light[1](Cur_light);

//Light l_list[2] = {light1, light2};
//Light l_list[1] = {Cur_light};
int l_num = 1;  //number of light sources you want to spawn 

vec4 MyWave1 = vec4(1.0, -2.0, 0.36, 0.3);
vec4 MyWave2 = vec4(0.0, 1.0, 0.12, 0.3);
vec4 MyWave3 = vec4(2.0, 1.0, 0.32, 0.3);
vec3 sphere_center_w1 = vec3(-0.7, -0.5, -2.5);
vec3 tangent_s = vec3(1, 0, 0);
vec3 binormal_s = vec3(0, 0, 1);
vec3 delta_center = vec3(0, 0, 0);
vec3 a = GerstnerWave(MyWave1, sphere_center_w1, tangent_s, binormal_s);
vec3 b = GerstnerWave(MyWave2, sphere_center_w1, tangent_s, binormal_s);
vec3 c = GerstnerWave(MyWave3, sphere_center_w1, tangent_s, binormal_s);
Sphere sphere_front = initSphere(vec3(0.1, 0.0, -1.5), 0.50, vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), 0.075, 0.075, 0.2, 20, 0.01, 0.8, 0.95);
Sphere sphere_back = initSphere(vec3(-0.7, -0.5, -2.5), 0.45, vec3(0.7, 0.7, 0.7), vec3(0.7, 0.7, 0.7), vec3(1.0, 1.0, 1.0), 0.15, 0.25, 1.0, 20, 0.75, 0.0, 0.0);
Sphere sphere_water1 = initSphere(sphere_center_w1 + a + b + c, 0.38, vec3(0.7, 0.7, 0.7), vec3(0.2, 0.9, 0.8), vec3(1.0, 1.0, 1.0), 0.15, 0.6, 0.35, 8, 0.0, 0.0, 0.0);
Sphere sphere_water2 = initSphere(vec3(0.5, -0.7, -2.0), 0.3, vec3(0.7, 0.7, 0.7), vec3(0.9, 0.2, 0.8), vec3(1.0, 1.0, 1.0), 0.15, 0.6, 0.35, 8, 0.0, 0.0, 0.0);
Sphere sphere_sun = initSphere(SunPos, 0.6, vec3(0.7, 0.7, 0.7), vec3(1.0, 0.3, 0.0), vec3(1.0, 1.0, 1.0), 0.15, 0.6, 0.35, 8, 0.0, 0.0, 0.0);
Sphere sphere_moon = initSphere(MoonPos, 0.6, vec3(0.7, 0.7, 0.7), vec3(0.4, 0.4, 0.1), vec3(1.0, 1.0, 1.0), 0.15, 0.6, 0.35, 8, 0.0, 0.0, 0.0);
//Sphere s_list[2] = {sphere_front, sphere_back};
Sphere s_list[2] = Sphere[2](sphere_water1, sphere_water2);
int s_num = 2;

float water_width = 40.0;
float water_depth = -0.4;
float floor_width = 40.0;
float floor_depth = -1.0;
Triangle Triangle_floor1 = initTriangle(vec3(-floor_width, floor_depth, -20.0), vec3(-floor_width, floor_depth, 0.0), vec3(floor_width, floor_depth, -20.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 0.8, 0.0), vec3(1.0, 1.0, 1.0), 0.15, 0.7, 0.25, 20, 0.0, 0.0, 0.0, false);
Triangle Triangle_floor2 = initTriangle(vec3(floor_width, floor_depth, 0.0), vec3(-floor_width, floor_depth, 0.0), vec3(floor_width, floor_depth, -20.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 0.8, 0.0), vec3(1.0, 1.0, 1.0), 0.15, 0.7, 0.25, 20, 0.0, 0.0, 0.0, false);
Triangle Triangle_water1 = initTriangle(vec3(-water_width, water_depth, -20.0), vec3(-water_width, water_depth, 10.0), vec3(water_width, water_depth, -20.0), vec3(1.0, 1.0, 1.0), vec3(0.2, 0.6, 0.4), vec3(1.0, 1.0, 1.0), 0.075, 0.05, 0.5, 40, 0.2, 0.7, 1.08, true);
Triangle Triangle_water2 = initTriangle(vec3(water_width, water_depth, 10.0), vec3(-water_width, water_depth, 10.0), vec3(water_width, water_depth, -20.0), vec3(1.0, 1.0, 1.0), vec3(0.2, 0.6, 0.4), vec3(1.0, 1.0, 1.0), 0.075, 0.05, 0.5, 40, 0.2, 0.7, 1.08, true);
Triangle t_list[4] = Triangle[4](Triangle_water1, Triangle_water2, Triangle_floor1, Triangle_floor2);
//Triangle t_list[2] = {Triangle_water1, Triangle_water2};
//Triangle t_list[2] = {Triangle_floor1, Triangle_floor2};
int t_num = 4;

/*Light sources and Primitive Shapes Init END*/

/*Local Illumination BEGIN*/
float SphereHit(Sphere sphere, Ray ray)
{
	vec3 oc = ray.origin - sphere.center;
	
	float a = dot(ray.direction, ray.direction);
	float b = 2.0 * dot(oc, ray.direction);
	float c = dot(oc, oc) - sphere.radius * sphere.radius;

	float discriminant = b * b - 4 * a * c;
    if(discriminant >= 0.0 && b <= 0){
        float t = (-b-sqrt(b*b-4*a*c))/(2*a);
        return t;
    }
    return 99999999.9;
}

float InternalSphereHit(Sphere sphere, Ray ray)
{
	vec3 oc = ray.origin - sphere.center;
	
	float a = dot(ray.direction, ray.direction);
	float b = 2.0 * dot(oc, ray.direction);
	float c = dot(oc, oc) - sphere.radius * sphere.radius;

	float discriminant = b * b - 4 * a * c;
    if(discriminant >= 0.0){
        float t = (-b+sqrt(b*b-4*a*c))/(2*a);
        return t;
    }
    return 99999999.9;
}

float TriangleHit(Triangle Triangle, Ray ray)
{
    vec3 dir = ray.direction;
    vec3 E1 = Triangle.v2 - Triangle.v1;
    vec3 E2 = Triangle.v3 - Triangle.v1;
    vec3 S = ray.origin - Triangle.v1;
    vec3 S1 = cross(dir, E2);
    vec3 S2 = cross(S, E1);
    
    if(dot(S1, E1) == 0){
        return 99999999.9;
    }
    else{
        float t = dot(S2, E2)/dot(S1, E1);
        float u = dot(S1, S)/dot(S1, E1);
        float v = dot(S2, dir)/dot(S1, E1);
        if(t > 0 && u > 0 && v > 0 && (1-v-u) > 0){
            return t;
        }
        else{
            return 99999999.9;
        }
    }
    
    //return true;
}

vec3 cal_SphereNormal(float t, Ray ray, Sphere sphere)
{
    vec3 hit = ExtendedRay(ray, t);
    vec3 s_normal = normalize(hit - sphere.center);
    return s_normal;
}

vec3 cal_wave_normal(vec3 normal){
    vec3 wave_normal;
    wave_normal = normal + normalize(vec3(1.0, 0.0, 0.0)) * sin(time_t);
    wave_normal = normalize(wave_normal);
    return wave_normal;
}


vec2 randomVec(vec2 uv)
{
	float tmp0 = dot(uv, vec2(127.1, 311.7));
    float tmp1 = dot(uv, vec2(327.2, 231.9));
	//float a = -1.0 + 2.0 * fract(sin(tmp0) * 43758.5453123);
    float a = fract(sin(tmp0) * 43758.5453123);
	float b = fract(sin(tmp1) * 85616.4864641);
    //return vec2(a, b);
    return vec2(a, b);
}

vec3 cal_circle_wave(vec3 gridPoint, vec3 centerPoint)
{
    vec2 circle_dir = vec2(vec3(gridPoint - centerPoint).x, vec3(gridPoint - centerPoint).z);
    float x = circle_dir.x;
    float z = circle_dir.y;
    float dis = sqrt(pow(x,2) + pow(z,2));
    vec3 dydx = vec3(1, -PAI*sin(PAI*(dis-time_t))*x/dis, 0);
    vec3 dydz = vec3(0, -PAI*sin(PAI*(dis-time_t))*z/dis, 1);
    if(dis < 0.001){
        dydx = vec3(1, PAI*cos(PAI*(dis-time_t)), 0);
        if(x < 0)
            dydx = vec3(1, -PAI*cos(PAI*(dis-time_t)), 0);
        dydz = vec3(0, PAI*cos(PAI*(dis-time_t)), 1);
        if(z < 0)
            dydz = vec3(0, -PAI*cos(PAI*(dis-time_t)), 1);
    }
    vec3 circle_normal = normalize(cross(dydx,dydz));
    if(circle_normal.y < 0)
        circle_normal = -circle_normal;
    return circle_normal;
}
float easecurve(float a, float b, float t)
{
    float k = 6.0 * pow(t, 5) - 15.0 * pow(t, 4) + 10.0 * pow(t, 3);
    float result = (1.0-k) * a + k * b;
    return result;
}

vec3 Gerstner_Circle_Wave (vec4 wave, vec3 p, inout vec3 tangent, inout vec3 binormal) 
{
    float steepness = wave.z * easecurve(1.0, 0.0, (time_t - timeStamp)/4);
    if((time_t - timeStamp)/4 > 1){
        steepness = 0;
    }
    float wavelength = wave.w;
    float k = 2 * PAI / wavelength;
    float c = sqrt(9.8 / k);                // _WaveSpeed
    vec2 dir = normalize(wave.xy);
    if(length(p.xz) > c * (time_t - timeStamp)/1.5){
        return vec3(0.0, 0.0, 0.0);
    }
    float f = k * (length(p.xz) - c * (time_t - timeStamp)/1.5);
    
    float a = steepness / k;                // _Amplitude

    tangent += vec3(-dir.x * dir.x * (steepness * sin(f)), dir.x * (steepness * cos(f)), -dir.x * dir.y * (steepness * sin(f)));
    binormal += vec3(-dir.x * dir.y * (steepness * sin(f)), dir.y * (steepness * cos(f)), -dir.y * dir.y * (steepness * sin(f)));

    return vec3(dir.x * (a * cos(f)), a * sin(f), dir.y * (a * cos(f)));
}

vec3 cal_TriangleNormal(Ray ray, inout Triangle triangle, float t)
{
    vec3 t_normal = cross((triangle.v1-triangle.v3), (triangle.v2-triangle.v3));
    t_normal = normalize(t_normal);
    if (dot(ray.direction, t_normal) > 0)
        t_normal = -t_normal;
    if(triangle.iswave){
        //GerstnerWave
        vec3 tangent = vec3(1, 0, 0);
        vec3 binormal = vec3(0, 0, 1);
        vec3 gridPoint = ExtendedRay(ray, t);
        vec4 NoiseWave = vec4(randomVec(gridPoint.xz), 0.2, 0.5);
        vec3 offset = vec3(0, 0, 0);
        offset += GerstnerWave(MyWave1, gridPoint, tangent, binormal);
        offset += GerstnerWave(MyWave2, gridPoint, tangent, binormal);
        offset += GerstnerWave(MyWave3, gridPoint, tangent, binormal);
        //GerstnerWave(NoiseWave, gridPoint, tangent, binormal);
        //t_normal = normalize(cross(binormal, tangent));
        
        //circle wave 1
        vec3 centerPoint = vec3(0.0, water_depth, -1.0);
        vec2 circle_dir = vec2(vec3(gridPoint - centerPoint).x, vec3(gridPoint - centerPoint).z);
        vec4 CircleWave = vec4(circle_dir, 0.3, 0.3);
        //Gerstner_Circle_Wave(CircleWave, gridPoint-centerPoint, tangent, binormal);
        vec3 centerPoint1 = WaterCenter;
        vec2 circle_dir1 = vec2(vec3(gridPoint - centerPoint1).x, vec3(gridPoint - centerPoint1).z);
        vec4 CircleWave1 = vec4(circle_dir1, 0.4, 0.4);
        offset += Gerstner_Circle_Wave(CircleWave1, gridPoint-centerPoint1, tangent, binormal);
        t_normal = normalize(cross(binormal, tangent));
        float maxsteep = MyWave1.z + MyWave2.z + MyWave3.z;
        triangle.kw = (maxsteep+offset.y)/(maxsteep-offset.y);
    }
        
    return t_normal;
}

vec3 cal_dir_lightsrc(float t, Ray ray, Light light)    //Light direction
{
    vec3 pTOlight = light.position -  ExtendedRay(ray, t);
    vec3 dir_lightsrc = normalize(pTOlight);
    return dir_lightsrc;
}

float cal_dis_lightsrc(float t, Ray ray, Light light)   //Light distance
{
    vec3 pTOlight = light.position -  ExtendedRay(ray, t);
    float dis_lightsrc = length(pTOlight);
    return dis_lightsrc;
}

vec3 cal_dir_rfl(vec3 dir_lightsrc, vec3 normal)     //Mirror reflect direction
{
    float cos_dn = dot(dir_lightsrc, normal);
    vec3 dir_rfl = cos_dn * normal * 2 - dir_lightsrc;
    dir_rfl = normalize(dir_rfl);
    return dir_rfl;
}

Ray cal_dir_trm_sph(float t, Ray ray, vec3 normal_in, Sphere sphere)
{
    Ray trm_in, trm_out;
    vec3 dir = ray.direction;
    float eta = sphere.eta;
    vec3 dir_trm_in, dir_trm_out;

    if(eta >= 1){
        dir_trm_in = (1/eta)* dir + ((1/eta)*(-dot(dir, normal_in))-sqrt(1+(pow((1/eta),2)*(pow(-dot(dir, normal_in),2)-1))))*normal_in;
        dir_trm_in = normalize(dir_trm_in);
        trm_in = initRay(ExtendedRay(ray, t), dir_trm_in);
        float t_in = InternalSphereHit(sphere, trm_in);
        vec3 normal_out = -cal_SphereNormal(t_in, trm_in, sphere);
        if(1-pow(dot(dir_trm_in, normal_out), 2) > pow((1/eta), 2)){
            dir_trm_out = cal_dir_rfl(-dir_trm_in, normal_out);
        }
        else{
            dir_trm_out = (eta)* dir_trm_in + ((eta)*(-dot(dir_trm_in, normal_out))-sqrt(1+(pow((eta),2)*(pow(-dot(dir_trm_in, normal_out),2)-1))))*normal_out;
            dir_trm_out = normalize(dir_trm_out);
        }
        trm_out = initRay(ExtendedRay(trm_in, t_in), dir_trm_out);
    }
    else{
        if(1-pow(dot(dir, normal_in), 2) > pow(eta, 2)){
            dir_trm_out = cal_dir_rfl(-dir, normal_in);
            trm_out = initRay(ExtendedRay(ray, t), dir_trm_out);
        }
        else{
            dir_trm_in = (1/eta)* dir + ((1/eta)*(-dot(dir, normal_in))-sqrt(1+(pow((1/eta),2)*(pow(-dot(dir, normal_in),2)-1))))*normal_in;
            trm_in = initRay(ExtendedRay(ray, t), dir_trm_in);
            float t_in = InternalSphereHit(sphere, trm_in);
            vec3 normal_out = -cal_SphereNormal(t_in, trm_in, sphere);
            dir_trm_out = (eta)* dir_trm_in + ((eta)*(-dot(dir_trm_in, normal_out))-sqrt(1+(pow((eta),2)*(pow(-dot(dir_trm_in, normal_out),2)-1))))*normal_out;
            trm_out = initRay(ExtendedRay(trm_in, t_in), dir_trm_out);
        }   
    }

    return trm_out;
}
Ray cal_dir_trm_tri(float t, Ray ray, vec3 normal_in, Triangle triangle)
{
    Ray trm_in, trm_out;
    vec3 dir = ray.direction;
    float eta = triangle.eta;
    vec3 dir_trm_in;

    if(eta >= 1){
        dir_trm_in = (1/eta)* dir + ((1/eta)*(-dot(dir, normal_in))-sqrt(1+(pow((1/eta),2)*(pow(-dot(dir, normal_in),2)-1))))*normal_in;
        dir_trm_in = normalize(dir_trm_in);
    }
    else{
        if(1-pow(dot(dir, normal_in), 2) > pow(eta, 2)){
            dir_trm_in = cal_dir_rfl(-dir, normal_in);
        }
        else{
            dir_trm_in = (1/eta)* dir + ((1/eta)*(-dot(dir, normal_in))-sqrt(1+(pow((1/eta),2)*(pow(-dot(dir, normal_in),2)-1))))*normal_in;
        }   
    }
    trm_in = initRay(ExtendedRay(ray, t) + 0.1 * dir_trm_in, dir_trm_in);

    return trm_in;
}

vec3 cal_Ambient(vec3 amb_c, float ka)
{
    vec3 Amb_Color = ka * amb_c;
    return Amb_Color;
}
vec3 cal_Diffuse(vec3 dif_c, vec3 normal, vec3 dir_lightsrc, float kd)
{
    vec3 Dif_Color = kd * dif_c * dot(normal, dir_lightsrc);
    return Dif_Color;
}
vec3 cal_Specular(Ray ray, vec3 dir_mirror, Light light, vec3 spec_c, float ks, int ke)
{
    vec3 Spec_Color = vec3(0.0, 0.0, 0.0);
    if(dot(-ray.direction, dir_mirror) > 0)
        Spec_Color = ks * spec_c * pow(dot(-ray.direction, dir_mirror), ke);
        //Spec_Color = specular;
        
    return Spec_Color;
}
bool is_shadow(float t_in, Ray ray, Light light)
{
    vec3 point = ExtendedRay(ray, t_in);
    vec3 dir = light.position- point;
    Ray p_ray = initRay(point, dir);
    p_ray.origin += 0.01 * p_ray.direction;
    int t_object[2] = int[2](0, 0);   //1 sphere, 2 triangle, {shape, num}
    float t_buffer = 9999999.9;
    float t = 9999999.9;
    if(t_num > 0){
        for(int i = 0; i < t_num; i ++){
            t = TriangleHit(t_list[i], p_ray);
            if(t < length(dir) && t < t_buffer){
                t_buffer = t;
                t_object[0] = 2;
                t_object[1] = i;
                //return true;
            } 
        }
    }
    if(s_num > 0){
        for(int i = 0; i < s_num; i ++){
            t = SphereHit(s_list[i], p_ray);
            if(t < length(dir) && t < t_buffer){
                t_buffer = t;
                t_object[0] = 1;
                t_object[1] = i;
                //return true;
            }        
        }
    }
    
    t = t_buffer;
    if (t_object[0] == 0){
        return false;
    }
    else if(t_object[0] == 1){
        if(s_list[t_object[1]].kt > 0.5){
            return false;
        }
        return true;
    }
    else if(t_object[0] == 2){
        if(t_list[t_object[1]].kt > 0.5){
            return false;
        }
        return true;
    }
    else{
        return true;
    }
}

/*Local Illumination END*/

/*Procedural Textures BEGIN*/

bool isodd(int n)
{
    int h = n/2;
    if(h*2 < n)
        return true;
    else
        return false;
}

vec3 concentric_circle_plane_xz(Triangle tri, Ray ray, float t)
{
    float r = 0.15;
    vec3 color = vec3(0.0, 0.0, 0.0);
    vec3 point = ExtendedRay(ray, t);
    float minrow = min(tri.v1.x, min(tri.v2.x, tri.v3.x));
    float maxrow = max(tri.v1.x, max(tri.v2.x, tri.v3.x));
    float mincol = min(tri.v1.z, min(tri.v2.z, tri.v3.z));
    float maxcol = max(tri.v1.z, max(tri.v2.z, tri.v3.z));
    vec3 centric = vec3((minrow+maxrow)/2, point.y, (mincol+maxcol)/2);
    float d = length(point - centric);
    int u = int(d/r);
    float angle = dot(normalize(point - centric), vec3(1.0, 0.0, 0.0));
    int v;
    if(angle > sqrt(3)/2 || (-sqrt(3)/2 < angle && angle < 0))
        v = 2;
    else 
        v = 1;
    if(isodd(u+v))
        color = vec3(1.0, 1.0, 0.0);
    else
        color = vec3(1.0, 0.0, 0.0);
    return color;
}


float MyPerlinNoise(vec2 uv)
{
    float noise = 0.0;
    float x0 = floor(uv.x);
    float x1 = ceil(uv.x);
    float y0 = floor(uv.y);
    float y1 = ceil(uv.y);
    float t0 = fract(uv.x);
    float t1 = fract(uv.y);
    float g00 = randomVec(vec2(x0, y0)).x;
    float g01 = randomVec(vec2(x0, y1)).x;
    float g10 = randomVec(vec2(x1, y0)).x;
    float g11 = randomVec(vec2(x1, y1)).x;
    float g1 = easecurve(g00, g10, t0);
    float g2 = easecurve(g01, g11, t0);
    noise = easecurve(g1, g2, t1);
    //float g1 = mix(g00, g10, t0);
    //float g2 = mix(g01, g11, t0);
    //noise = mix(g1, g2, t1);
    return noise;
}

vec3 perlin_plane_xz(Triangle tri, Ray ray, float t)
{
    vec3 color = vec3(0.0, 0.0, 0.0);
    vec3 point = ExtendedRay(ray, t);
    float minrow = min(tri.v1.x, min(tri.v2.x, tri.v3.x));
    float maxrow = max(tri.v1.x, max(tri.v2.x, tri.v3.x));
    float mincol = min(tri.v1.z, min(tri.v2.z, tri.v3.z));
    float maxcol = max(tri.v1.z, max(tri.v2.z, tri.v3.z));
    float l_row = maxrow - minrow;
    float l_col = maxcol - mincol;
    float row_a = 0.25 * 1;                 //pattern row length
    float col_a = float(1.0/6.0) * 1;      //pattern col length
    //float u = (point.x-minrow)/l_row;
    //float v = (point.z-mincol)/l_col; 
    float u = (point.x-minrow)/row_a;
    float v = (point.z-mincol)/col_a;  
    float noise = MyPerlinNoise(vec2(u, v));
    color = vec3(0.3+noise*0.5, 0.24+noise*0.4, 0.0);
    return color;
}

vec3 checkerboard_plane_xz(Triangle tri, Ray ray, float t)
{
    float row_a = 0.25;                 //pattern row length
    float col_a = float(1.0/6.0);      //pattern col length
    vec3 color = vec3(0.0, 0.0, 0.0);
    vec3 point = ExtendedRay(ray, t);
    float minrow = min(tri.v1.x, min(tri.v2.x, tri.v3.x));
    float maxrow = max(tri.v1.x, max(tri.v2.x, tri.v3.x));
    float mincol = min(tri.v1.z, min(tri.v2.z, tri.v3.z));
    float maxcol = max(tri.v1.z, max(tri.v2.z, tri.v3.z));
    float u = (point.x-minrow)/row_a;
    float v = (point.z-mincol)/col_a;  
    int sum = int(floor(u) + floor(v));
    if(isodd(sum))
        color = vec3(1.0, 1.0, 0.0);
    else
        color = vec3(1.0, 0.0, 0.0);
    return color;
}

/*Procedural Textures END*/

//Caustics calculation BEGIN
float cal_caustics(float t, Ray ray, vec3 normal, vec3 dir_lightsrc)
{
    //two parameters to tweak the caustics
    float area_coeff = 2.5;//the higher, the more "blurry" the caustics gets, the less area has caustics
    float project_length = 2;//the higher, the higher contrast in area, the higher the diffuse result

    vec3 hit_point = ray.origin + t * ray.direction;
    vec3 tris_below[3];
    tris_below[0] = hit_point + cross(normal, vec3(1,0,0)) * 2 * sqrt(3.000) / 3 * area_coeff;
    tris_below[1] = hit_point + (-cross(normal, vec3(1,0,0)) * sqrt(3.000) / 3 - vec3(1,0,0)) * area_coeff;
    tris_below[2] = hit_point + (-cross(normal, vec3(1,0,0)) * sqrt(3.000) / 3 + vec3(1,0,0)) * area_coeff;
    float area_below = sqrt(3.000) * area_coeff;


    vec3 tris_above[3];
    for(int i = 0; i < 3; i++) 
    {
        Ray caus_ray = initRay(tris_below[i], dir_lightsrc);
        
        //ray - triangle intesection, get the refraction position and which triangle
        float dist = 99999999.9;
        int tri_index = -1;
        for(int j = 0; j < t_list.length(); j++)
        {
            if(!t_list[j].iswave) continue;

            if( dist > TriangleHit(t_list[j], caus_ray)) 
            {
                dist = TriangleHit(t_list[j], caus_ray);
                tri_index = j;
            }
        }

        //if ray intersect with no wave triangle, which happens when the sun/moon is very low on the horizon
        if(tri_index == -1) return 1;
        
        vec3 normal_in = cal_TriangleNormal(caus_ray, t_list[tri_index], dist);
        caus_ray = cal_dir_trm_tri(dist, caus_ray, normal_in, t_list[tri_index]);
        tris_above[i] = caus_ray.origin + caus_ray.direction * project_length;//get the positions of triangles above water
    }
    vec3 AB = tris_above[1] - tris_above[0];
    vec3 AC = tris_above[2] - tris_above[0];
    float area_above = 0.5 * length(cross(AB, AC));

    return area_above / area_below /2.0;
}
//Caustics calculation END

/*Global Illumination BEGIN*/
vec3 Local_Illuminate(Ray ray)
{
    vec3 Color = Skycolor;
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = int[2](0, 0);   //1 sphere, 2 triangle, {shape, num}
    if(s_num > 0){
        for(int i = 0; i < s_num; i ++){    //multiple sphere
            t = SphereHit(s_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 1;
                t_object[1] = i;
            }
        }
    }
    if(t_num > 0){
        for(int i = 0; i < t_num; i ++){    //multiple triangle
            t = TriangleHit(t_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 2;
                t_object[1] = i;
            }
        }
    }
    
    int sid = t_object[1];  //shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    t = tBuffer;   
    switch (t_object[0]){
        case 0:
            Color = Skycolor;
            break;
        case 1: // Sphere
            Color = vec3(0.0, 0.0, 0.0);
            cur_amb = s_list[sid].amb_c;
            cur_dif = s_list[sid].dif_c;
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_SphereNormal(t, ray, s_list[sid]);
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, s_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);
                    
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd * kc) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks * kc, s_list[sid].ke)) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            //vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            //vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            cur_amb = t_list[sid].amb_c;
            //cur_amb = tex1;
            if(!t_list[sid].iswave)
                cur_dif = tex3;
            else
                cur_dif = t_list[sid].dif_c;
            cur_spec = t_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid], t); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);

                    Color += (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd * kc)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks * kc, t_list[sid].ke)) / l_num;
                } 
            }
            break;
        default:
            break;
    }
     
    return Color;
}
vec3 RayTrace2(Ray ray)
{
    vec3 Color = Skycolor;
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = int[2](0, 0);   //1 sphere, 2 triangle, {shape, num}
    if(s_num > 0){
        for(int i = 0; i < s_num; i ++){    //multiple sphere
            t = SphereHit(s_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 1;
                t_object[1] = i;
            }
        }
    }
    if(t_num > 0){
        for(int i = 0; i < t_num; i ++){    //multiple triangle
            t = TriangleHit(t_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 2;
                t_object[1] = i;
            }
        }
    }
    
    int sid = t_object[1];  //shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    int MAX_DEPTH = 2;
    t = tBuffer;   
    switch (t_object[0]){
        case 0:
            Color = Skycolor;
            break;
        case 1: // Sphere
            Color = vec3(0.0, 0.0, 0.0);
            cur_amb = s_list[sid].amb_c;
            cur_dif = s_list[sid].dif_c;
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_SphereNormal(t, ray, s_list[sid]);
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, s_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);
                    
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd * kc) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks * kc, s_list[sid].ke)) / l_num;
                }
                if(s_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(-ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    Color += s_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    //Color += s_list[sid].kr * RayTrace(rfl) / l_num;
                }
                if(s_list[sid].kt > 0){
                    Ray trm = cal_dir_trm_sph(t, ray, normal, s_list[sid]);
                    Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    //Color += s_list[sid].kt * RayTrace(trm) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            //vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            //vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            cur_amb = t_list[sid].amb_c;
            //cur_amb = tex1;
            if(!t_list[sid].iswave)
                cur_dif = tex3;
            else
                cur_dif = t_list[sid].dif_c;
            cur_spec = t_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid], t); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);
                    
                    Color += (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd * kc)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks * kc, t_list[sid].ke)) / l_num;
                } 
                if(t_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t-0.01), dir_rfl);
                    Color += t_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    //Color += t_list[sid].kr * RayTrace(rfl) / l_num;
                }
                if(t_list[sid].kt > 0){
                    Ray trm = cal_dir_trm_tri(t, ray, normal, t_list[sid]);
                    //Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    Color += t_list[sid].kt * Local_Illuminate(trm) / l_num;
                }
            }
            break;
        default:
            break;
    }
     
    return Color;
}
vec3 RayTrace1(Ray ray)
{
    vec3 Color = Skycolor;
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = int[2](0, 0);   //1 sphere, 2 triangle, {shape, num}
    if(s_num > 0){
        for(int i = 0; i < s_num; i ++){    //multiple sphere
            t = SphereHit(s_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 1;
                t_object[1] = i;
            }
        }
    }
    if(t_num > 0){
        for(int i = 0; i < t_num; i ++){    //multiple triangle
            t = TriangleHit(t_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 2;
                t_object[1] = i;
            }
        }
    }
     
    int sid = t_object[1];  //shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    int MAX_DEPTH = 2;
    t = tBuffer;   
    switch (t_object[0]){
        case 0:
            Color = Skycolor;
            break;
        case 1: // Sphere
            Color = vec3(0.0, 0.0, 0.0);
            cur_amb = s_list[sid].amb_c;
            cur_dif = s_list[sid].dif_c;
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_SphereNormal(t, ray, s_list[sid]);
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, s_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);
                    
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd * kc) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks * kc, s_list[sid].ke)) / l_num;
                }
                if(s_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(-ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    //Color += s_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    Color += s_list[sid].kr * RayTrace2(rfl) / l_num;
                }
                if(s_list[sid].kt > 0){
                    Ray trm = cal_dir_trm_sph(t, ray, normal, s_list[sid]);
                    //Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    Color += s_list[sid].kt * RayTrace2(trm) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            //vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            //vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            cur_amb = t_list[sid].amb_c;
            //cur_amb = tex1;
            if(!t_list[sid].iswave)
                cur_dif = tex3;
            else
                cur_dif = t_list[sid].dif_c;
            cur_spec = t_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid], t); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);
                    
                    Color += (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd * kc)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks * kc, t_list[sid].ke)) / l_num;
                } 
                if(t_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t-0.01), dir_rfl);
                    //Color += t_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    Color += t_list[sid].kr * RayTrace2(rfl) / l_num;
                }
                if(t_list[sid].kt > 0){
                    Ray trm = cal_dir_trm_tri(t, ray, normal, t_list[sid]);
                    //Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    Color += t_list[sid].kt * RayTrace2(trm) / l_num;
                }
            }
            break;
        default:
            break;
    }
     
    return Color;
}
vec3 RayTrace(Ray ray)
{
    vec3 Color = Skycolor;
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = int[2](0, 0);   //1 sphere, 2 triangle, {shape, num}
    if(s_num > 0){
        for(int i = 0; i < s_num; i ++){    //multiple sphere
            t = SphereHit(s_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 1;
                t_object[1] = i;
            }
        }
    }
    if(t_num > 0){
        for(int i = 0; i < t_num; i ++){    //multiple triangle
            t = TriangleHit(t_list[i], ray);
            if(t < tBuffer){
                tBuffer = t;
                t_object[0] = 2;
                t_object[1] = i;
            }
        }
    }
    
    int sid = t_object[1];  //shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    int MAX_DEPTH = 2;
    t = tBuffer;   
    switch (t_object[0]){
        case 0:
            Color = Skycolor;
            break;
        case 1: // Sphere
            Color = vec3(0.0, 0.0, 0.0);
            cur_amb = s_list[sid].amb_c;
            cur_dif = s_list[sid].dif_c;
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_SphereNormal(t, ray, s_list[sid]);
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, s_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);
                    
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd * kc) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks * kc, s_list[sid].ke)) / l_num;
                }
                if(s_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(-ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    //Color += s_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    Color += s_list[sid].kr * RayTrace1(rfl) / l_num;
                }
                if(s_list[sid].kt > 0){
                    Ray trm = cal_dir_trm_sph(t, ray, normal, s_list[sid]);
                    //Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    Color += s_list[sid].kt * RayTrace1(trm) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            //vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            //vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            cur_amb = t_list[sid].amb_c;
            //cur_amb = tex1;
            if(!t_list[sid].iswave)
                cur_dif = tex3;
            else
                cur_dif = t_list[sid].dif_c;
            cur_spec = t_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid], t); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += t_list[sid].kw * cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    //caustics
                    float kc = 1.00; //caustics coefficient
                    if(!t_list[sid].iswave) kc = cal_caustics(t, ray, normal, dir_lightsrc);
                    
                    Color += t_list[sid].kw * (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd * kc)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks * kc, t_list[sid].ke)) / l_num;
                } 
                if(t_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    Color += t_list[sid].kw * t_list[sid].kr * Skycolor / l_num;
                    //Color += t_list[sid].kr * RayTrace1(rfl) / l_num;
                }
                if(t_list[sid].kt > 0){
                    Ray trm = cal_dir_trm_tri(t, ray, normal, t_list[sid]);
                    Color += t_list[sid].kt * Local_Illuminate(trm) / l_num;
                    //Color += t_list[sid].kt * RayTrace1(trm) / l_num;
                }
            }
            break;
        default:
            break;
    }
     
    return Color;
}

vec3 draw_light(Ray ray, Sphere light)
{
    vec3 Color = Skycolor;
    float t = SphereHit(light, ray);
    if(t < 99999999.9){
        Color = light.dif_c;
    }
    return Color;
}
/*Global Illumination END*/

void main()
{
    float u = screenCoord.x;
	float v = screenCoord.y;
	
	vec3 origin = vec3(0.0, 0.0, 0.0);
	vec3 direction1 = origin + vec3(-1.0, -1.0, -1.0) + u * vec3(2.0, 0.0, 0.0) + v * vec3(0.0, 2.0, 0.0);
    
    //Super-sampling
    Ray ray1 = initRay(origin, direction1);
    vec3 direction2 = origin + vec3(-1.0, -1.0, -1.0) + (u+0.0005) * vec3(2.0, 0.0, 0.0) + v * vec3(0.0, 2.0, 0.0);
    Ray ray2 = initRay(origin, direction2);
    vec3 direction3 = origin + vec3(-1.0, -1.0, -1.0) + u * vec3(2.0, 0.0, 0.0) + (v+0.0005) * vec3(0.0, 2.0, 0.0);
    Ray ray3 = initRay(origin, direction3);
    vec3 direction4 = origin + vec3(-1.0, -1.0, -1.0) + (u+0.0005) * vec3(2.0, 0.0, 0.0) + (v+0.0005) * vec3(0.0, 2.0, 0.0);
    Ray ray4 = initRay(origin, direction4);
    vec3 WaterDirecrion = origin + vec3(-1.0, -1.0, -1.0) + mousePos.x * vec3(2.0, 0.0, 0.0) + mousePos.y * vec3(0.0, 2.0, 0.0);
    Ray waterRay = initRay(origin, WaterDirecrion);
    Ray r_list[4] = Ray[4](ray1, ray2, ray3, ray4);
    int r_num = 1;  //number of ray you want to spawn 

    bool isSun = false;
    bool isMoon = false;
    vec3 FinalColor = Skycolor;
    float water_t1 = 99999.9;
    float water_t2 = 99999.9;
    water_t1 = TriangleHit(Triangle_water1, waterRay);
    water_t2 = TriangleHit(Triangle_water2, waterRay);
    if(water_t1 < 99999){
        WaterCenter = ExtendedRay(waterRay, water_t1);
    }
    else if(water_t2 < 99999){
        WaterCenter = ExtendedRay(waterRay, water_t2);
    }
    FinalColor = vec3(0.0, 0.0, 0.0);
    for(int i = 0; i < r_num; i ++){ 
        if(draw_light(r_list[i], sphere_sun) != Skycolor && SunPos.y > -0.4){
            FinalColor += draw_light(r_list[i], sphere_sun)/r_num;
            isSun = true;
        }  
        if(draw_light(r_list[i], sphere_moon) != Skycolor && MoonPos.y > -0.4){
            FinalColor += draw_light(r_list[i], sphere_moon)/r_num;
            isMoon = true;
        }      
        if(!isSun && !isMoon){
            if(SunPos.y > 0)
                FinalColor += (0.25+0.75*sin(DayTime))*RayTrace(r_list[i])/r_num; 
            else
                FinalColor += 0.25*RayTrace(r_list[i])/r_num; 
        }        
    }
    //FinalColor += RayTrace(ray1)/r_num;
    //FinalColor = cal_circle_wave(vec3(-1+2*u, -1, -1+2*v), vec3(0, -1, 0));
    FragColor = vec4(FinalColor, 1.0);
    //Showing click on screen

    // float clickPersision = 0.005f;
    // if(length(mousePos.xy - screenCoord) < clickPersision){
    //     FragColor = vec4(1,1,1,1);
    // }
}