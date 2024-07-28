#version 330 core
 
in vec3 ourColor;
in vec2 screenCoord;
out vec4 FragColor;

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

Triangle initTriangle(vec3 v1, vec3 v2, vec3 v3, vec3 amb_c, vec3 dif_c, vec3 spec_c, float ka, float kd, float ks, int ke, float kr, float kt, float eta)
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
    return triangle;
}

/*Light sources and Primitive Shapes Init BEGIN*/

Light light1 = initLight(vec3(0.0, 5.0, 1.0), vec3(1.0, 1.0, 1.0));
Light light2 = initLight(vec3(2.0, 5.0, -1.0), vec3(1.0, 1.0, 1.0));
Light l_list[2] = {light1, light2};
int l_num = 1;  //number of light sources you want to spawn 

Sphere sphere_front = initSphere(vec3(0.1, 0.0, -1.5), 0.50, vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), 0.075, 0.075, 0.2, 20, 0.01, 0.8, 0.95);
Sphere sphere_back = initSphere(vec3(-0.7, -0.5, -2.5), 0.45, vec3(0.7, 0.7, 0.7), vec3(0.7, 0.7, 0.7), vec3(1.0, 1.0, 1.0), 0.15, 0.25, 1.0, 20, 0.75, 0.0, 0.0);
Sphere s_list[2] = {sphere_front, sphere_back};
int s_num = 2;

Triangle Triangle_floor1 = initTriangle(vec3(-2.1, -1.0, -1.0), vec3(-2.1, -1.0, -4.5), vec3(0.9, -1.0, -1.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 1.0, 1.0), 0.15, 0.6, 0.35, 8, 0.0, 0.0, 0.0);
Triangle Triangle_floor2 = initTriangle(vec3(0.9, -1.0, -4.5), vec3(-2.1, -1.0, -4.5), vec3(0.9, -1.0, -1.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 0.0), vec3(1.0, 1.0, 1.0), 0.15, 0.6, 0.35, 8, 0.0, 0.0, 0.0);
Triangle t_list[2] = {Triangle_floor1, Triangle_floor2};
int t_num = 2;

/*Light sources and Primitive Shapes Init END*/

/*Local and Global Illumination related functions BEGIN*/
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

vec3 cal_TriangleNormal(Ray ray, Triangle triangle)
{
    vec3 t_normal = cross((triangle.v1-triangle.v3), (triangle.v2-triangle.v3));
    t_normal = normalize(t_normal);
    if (dot(ray.direction, t_normal) > 0)
        t_normal = -t_normal;
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

vec3 cal_dir_rfl(vec3 dir_lightsrc, vec3 normal)     //Calculate reflection direction of all shapes
{
    float cos_dn = dot(dir_lightsrc, normal);
    vec3 dir_rfl = cos_dn * normal * 2 - dir_lightsrc;
    dir_rfl = normalize(dir_rfl);
    return dir_rfl;
}

Ray cal_dir_trm(float t, Ray ray, vec3 normal_in, Sphere sphere)    //calculate sphere transmission out ray direction
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
    int t_object[2] = {0, 0};   //1 sphere, 2 triangle, {shape, num}
    float t_buffer = 9999999.9;
    float t = 9999999.9;
    for(int i = 0; i < t_num; i ++){
        t = TriangleHit(t_list[i], p_ray);
        if(t < length(dir) && t < t_buffer){
            t_buffer = t;
            t_object[0] = 2;
            t_object[1] = i;
            //return true;
        } 
    }
    for(int i = 0; i < s_num; i ++){
        t = SphereHit(s_list[i], p_ray);
        if(t < length(dir) && t < t_buffer){
            t_buffer = t;
            t_object[0] = 1;
            t_object[1] = i;
            //return true;
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
    else{
        return true;
    }
}

/*Local and Global Illumination related functions END*/

/*Procedural Textures(Checkpoint 4) BEGIN*/ 

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

float randomVec(vec2 uv)
{
	float tmp0 = dot(uv, vec2(127.1, 311.7));
    float tmp1 = dot(uv, vec2(327.2, 231.9));
	//float a = -1.0 + 2.0 * fract(sin(tmp0) * 43758.5453123);
    float a = fract(sin(tmp0) * 43758.5453123);
	float b = fract(sin(tmp1) * 85616.4864641);
    //return vec2(a, b);
    return a;
}

float easecurve(float a, float b, float t)
{
    float k = 6.0 * pow(t, 5) - 15.0 * pow(t, 4) + 10.0 * pow(t, 3);
    float result = (1.0-k) * a + k * b;
    return result;
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
    float g00 = randomVec(vec2(x0, y0));
    float g01 = randomVec(vec2(x0, y1));
    float g10 = randomVec(vec2(x1, y0));
    float g11 = randomVec(vec2(x1, y1));
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
    float row_a = 0.25;                 //pattern row length
    float col_a = float(1.0/6.0);      //pattern col length
    //float u = (point.x-minrow)/l_row;
    //float v = (point.z-mincol)/l_col; 
    float u = (point.x-minrow)/row_a;
    float v = (point.z-mincol)/col_a;  
    float noise = MyPerlinNoise(vec2(u, v));
    color = vec3(0.0, noise, 1.0);
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

/*Procedural Textures(Checkpoint 4) END*/

/*Local & Global Illumination BEGIN*/
vec3 Local_Illuminate(Ray ray)  // Fourth(and last) ray detection path
{
    vec3 Color = vec3(0.1, 0.6, 0.9);
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = {0, 0};   // {shape, index} shape: 0 no hit, 1 sphere, 2 triangle index:the index number of the shape hit in the list
    for(int i = 0; i < s_num; i ++){    //hit detect through all the spheres
        t = SphereHit(s_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 1;
            t_object[1] = i;
        }
    }
    for(int i = 0; i < t_num; i ++){    //hit detect through all the triangles
        t = TriangleHit(t_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 2;
            t_object[1] = i;
        }
    }
    
    int sid = t_object[1];  //get closest hit shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    t = tBuffer;   //get closest hit distance
    switch (t_object[0]){
        case 0:
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
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks, s_list[sid].ke)) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            //choose the texture that you want to use on the plane
            //cur_amb = s_list[sid].amb_c;
            cur_amb = tex1;
            //cur_dif = s_list[sid].dif_c;
            cur_dif = tex1;
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid]); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    Color += (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks, t_list[sid].ke)) / l_num;
                } 
            }
            break;
        default:
            break;
    }
     
    return Color;
}
vec3 RayTrace2(Ray ray)     // Third ray detection path
{
    vec3 Color = vec3(0.1, 0.6, 0.9);
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = {0, 0};   // {shape, index} shape: 0 no hit, 1 sphere, 2 triangle index:the index number of the shape hit in the list
    for(int i = 0; i < s_num; i ++){    //hit detect through all the spheres
        t = SphereHit(s_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 1;
            t_object[1] = i;
        }
    }
    for(int i = 0; i < t_num; i ++){    //hit detect through all the triangles
        t = TriangleHit(t_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 2;
            t_object[1] = i;
        }
    }
    
    int sid = t_object[1];  //get closest hit shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    int MAX_DEPTH = 2;
    t = tBuffer;   //get closest hit distance
    switch (t_object[0]){
        case 0:
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
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks, s_list[sid].ke)) / l_num;
                }
                if(s_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(-ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    Color += s_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    //Color += s_list[sid].kr * RayTrace(rfl) / l_num;
                }
                if(s_list[sid].kt > 0){
                    Ray trm = cal_dir_trm(t, ray, normal, s_list[sid]);
                    Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    //Color += s_list[sid].kt * RayTrace(trm) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            //choose the texture that you want to use on the plane
            //cur_amb = s_list[sid].amb_c;
            cur_amb = tex1;
            //cur_dif = s_list[sid].dif_c;
            cur_dif = tex1;
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid]); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    Color += (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks, t_list[sid].ke)) / l_num;
                } 
                if(t_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    Color += t_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    //Color += t_list[sid].kr * RayTrace(rfl) / l_num;
                }
                if(t_list[sid].kt > 0){
                    //Ray trm = cal_dir_trm(t, ray, normal, t_list[sid]);
                    //Color += t_list[sid].kt * Local_Illuminate(trm) / l_num;
                    //Color += t_list[sid].kt * RayTrace(trm) / l_num;
                }
            }
            break;
        default:
            break;
    }
     
    return Color;
}
vec3 RayTrace1(Ray ray) // Second ray detection path
{
    vec3 Color = vec3(0.1, 0.6, 0.9);
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = {0, 0};   // {shape, index} shape: 0 no hit, 1 sphere, 2 triangle index:the index number of the shape hit in the list
    for(int i = 0; i < s_num; i ++){    //hit detect through all the spheres
        t = SphereHit(s_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 1;
            t_object[1] = i;
        }
    }
    for(int i = 0; i < t_num; i ++){    //hit detect through all the triangles
        t = TriangleHit(t_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 2;
            t_object[1] = i;
        }
    }
    
    int sid = t_object[1];  //get closest hit shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    int MAX_DEPTH = 2;
    t = tBuffer;   //get closest hit distance
    switch (t_object[0]){
        case 0:
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
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks, s_list[sid].ke)) / l_num;
                }
                if(s_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(-ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    //Color += s_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    Color += s_list[sid].kr * RayTrace2(rfl) / l_num;
                }
                if(s_list[sid].kt > 0){
                    Ray trm = cal_dir_trm(t, ray, normal, s_list[sid]);
                    //Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    Color += s_list[sid].kt * RayTrace2(trm) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            //choose the texture that you want to use on the plane

            //cur_amb = s_list[sid].amb_c;
            cur_amb = tex1;
            //cur_dif = s_list[sid].dif_c;
            cur_dif = tex1;
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid]); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    Color += (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks, t_list[sid].ke)) / l_num;
                } 
                if(t_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    //Color += t_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    Color += t_list[sid].kr * RayTrace2(rfl) / l_num;
                }
                if(t_list[sid].kt > 0){
                    //Ray trm = cal_dir_trm(t, ray, normal, t_list[sid]);
                    //Color += t_list[sid].kt * Local_Illuminate(trm) / l_num;
                    //Color += t_list[sid].kt * RayTrace2(trm) / l_num;
                }
            }
            break;
        default:
            break;
    }
     
    return Color;
}
vec3 RayTrace(Ray ray)  // First ray detection path
{
    vec3 Color = vec3(0.1, 0.6, 0.9);
    float t = 0;
    float tBuffer = 99999999.9;
    int t_object[2] = {0, 0};   // {shape, index} shape: 0 no hit, 1 sphere, 2 triangle index:the index number of the shape hit in the list
    for(int i = 0; i < s_num; i ++){    //hit detect through all the spheres
        t = SphereHit(s_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 1;
            t_object[1] = i;
        }
    }
    for(int i = 0; i < t_num; i ++){    //hit detect through all the triangles
        t = TriangleHit(t_list[i], ray);
        if(t < tBuffer){
            tBuffer = t;
            t_object[0] = 2;
            t_object[1] = i;
        }
    }
    
    int sid = t_object[1];  //get closest hit shape index
    vec3 cur_amb = vec3(0.0, 0.0, 0.0);
    vec3 cur_dif = vec3(0.0, 0.0, 0.0);
    vec3 cur_spec = vec3(0.0, 0.0, 0.0);
    int MAX_DEPTH = 2;
    t = tBuffer;   //get closest hit distance
    switch (t_object[0]){
        case 0:
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
                    Color += (cal_Ambient(cur_amb, s_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, s_list[sid].kd) 
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, s_list[sid].ks, s_list[sid].ke)) / l_num;
                }
                if(s_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(-ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    //Color += s_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    Color += s_list[sid].kr * RayTrace1(rfl) / l_num;
                }
                if(s_list[sid].kt > 0){
                    Ray trm = cal_dir_trm(t, ray, normal, s_list[sid]);
                    //Color += s_list[sid].kt * Local_Illuminate(trm) / l_num;
                    Color += s_list[sid].kt * RayTrace1(trm) / l_num;
                }
            }
            break;
        case 2: //Triangle
            Color = vec3(0.0, 0.0, 0.0);
            vec3 tex1 = checkerboard_plane_xz(t_list[sid], ray, t);
            vec3 tex2 = concentric_circle_plane_xz(t_list[sid], ray, t);
            vec3 tex3 = perlin_plane_xz(t_list[sid], ray, t);
            //cur_amb = s_list[sid].amb_c;
            cur_amb = tex1;     //choose the texture that you want to use on the plane
            //cur_dif = s_list[sid].dif_c;
            cur_dif = tex1;     //choose the texture that you want to use on the plane
            cur_spec = s_list[sid].spec_c;
            for(int j = 0; j < l_num; j ++){    //multiple light
                vec3 normal = cal_TriangleNormal(ray, t_list[sid]); 
                vec3 dir_lightsrc = cal_dir_lightsrc(t, ray, l_list[j]);
                vec3 dir_mirror = cal_dir_rfl(dir_lightsrc, normal);
                if(is_shadow(t, ray, l_list[j])){
                    Color += cal_Ambient(cur_amb, t_list[sid].ka) / l_num;
                }
                else{
                    Color += (cal_Ambient(cur_amb, t_list[sid].ka) 
                          + cal_Diffuse(cur_dif, normal, dir_lightsrc, t_list[sid].kd)
                          + cal_Specular(ray, dir_mirror, l_list[j], cur_spec, t_list[sid].ks, t_list[sid].ke)) / l_num;
                } 
                if(t_list[sid].kr > 0){
                    vec3 dir_rfl = cal_dir_rfl(ray.direction, normal);
                    Ray rfl = initRay(ExtendedRay(ray, t), dir_rfl);
                    //Color += t_list[sid].kr * Local_Illuminate(rfl) / l_num;
                    Color += t_list[sid].kr * RayTrace1(rfl) / l_num;
                }
                if(t_list[sid].kt > 0){
                    //Ray trm = cal_dir_trm(t, ray, normal, t_list[sid]);
                    //Color += t_list[sid].kt * Local_Illuminate(trm) / l_num;
                    //Color += t_list[sid].kt * RayTrace1(trm) / l_num;
                }
            }
            break;
        default:
            break;
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

    Ray r_list[4] = {ray1, ray2, ray3, ray4};
    int r_num = 1;  //number of ray you want to spawn 

    for(int i = 0; i < r_num; i ++){
        FragColor.xyz += RayTrace(r_list[i])/r_num;  
        //FragColor.xyz += Local_Illuminate(r_list[i])/r_num;  
    }
    
	FragColor.w = 1.0;
    
}