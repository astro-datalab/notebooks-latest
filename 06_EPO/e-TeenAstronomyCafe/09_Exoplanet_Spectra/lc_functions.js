const data = source.data;
const rad = r.value; // planet ardius
const x = data['x']
const y = data['y']

const data_p = source_planet.data;
const r_show = data_p['r'];
const x_p = data_p['x']
const y_p = data_p['y']

const t_now = t.value;
const t_show = data_p['time_now']
const f_show = data_p['flux_now']
const a_o_r = 6.0 // semi-major axis
const b = b_imp.value // impact parameter
const period = 24. // orbital period
const u_lin = 0.2 // linear limb darkening

var x_proj;
var y_proj;
var z;
var r_p_rs;
var f;
var theta1;
var theta2;
var z_fun;
var r_fun;
var mu;
var C;
var x1;
var Aint;

function area_intersect(z_fun,r_fun) {
      if (z_fun >= 1.0 + r_fun) {
            f = 0.0
      } else if (z_fun < 1.0 - r_fun) {
            f = r_fun * r_fun
      } else {
            x1 = (1.0 - r_fun * r_fun + z_fun * z_fun) / (2.0 * z_fun)
            theta1 = Math.acos(x1)
            theta2 = Math.acos((z_fun - x1)/r_fun)
            Aint = theta1 + theta2 * r_fun * r_fun - Math.sqrt(1.0 - x1 * x1) * z_fun
            f = Aint / Math.PI
      }
      return f
}

// approximate limb darkening by ignoring variations over planet
function limb_dark(z_fun,r_fun) {
      C = 1.0 / (1.0 - u_lin / 6.0)
      if (z_fun >= 1.0 + r_fun) {
            f = 1.0
      } else if (z_fun < 1.0) {
            mu = Math.sqrt(1.0 - z_fun * z_fun)
            f = C * (1.0 - u_lin * (1.0 - mu))
      } else {
            f = C * (1.0 - u_lin)
      }
      return f
}

function lcFunction(t_fun, rad_fun) {
      x_proj = a_o_r * Math.sin(t_fun * 2.0 * Math.PI / period)
      y_proj = b * Math.cos(t_fun * 2.0 * Math.PI / period)
      z = Math.sqrt(x_proj * x_proj + y_proj * y_proj)
      r_p_rs = rad_fun / 10.0
      f = 1.0 - area_intersect(z, r_p_rs) * limb_dark(z,r_p_rs)
      return f * 100.0
}

r_show[0] = rad * 1.0;
x_p[0] = Math.sin(t_now * 6.28 / period) * a_o_r * 10.0;
y_p[0] = Math.cos(t_now * 6.28 / period) * b * 10.0;
t_show[0] = t_now;
f_show[0] = lcFunction(t_now, rad);

for (var i = 0; i < x.length; i++) {
    y[i] = lcFunction(x[i], rad);
}

source.change.emit();
source_planet.change.emit();
