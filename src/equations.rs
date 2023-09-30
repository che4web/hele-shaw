use pdifflib::field::{Field};
use pdifflib::finit_diff::poisson_relax;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use crate::PEREODIC;

#[derive(Serialize, Deserialize, Copy, Clone,Debug)]
pub struct ParamsConc {
    pub le: f64,
    pub sor:f64,
    #[serde(default ="zero")]
    pub sor1:f64,
    pub bm:f64,
    pub l_sed:f64,
}
fn zero() ->f64{
    0.0
}
trait Equation{
    fn write_map(&mut self, time:f64) ;

}

#[derive(Serialize, Deserialize, Copy, Clone,Debug)]
pub struct Params {
    pub rel: f64,
    pub gr_v: f64,
    pub omega: f64,
    pub pr: f64,
    pub time: f64,
    pub params_c: ParamsConc,
}
use crate::{H, NX, NY};
pub struct Temperatura {
    pub f: Array2<f64>,
    delta: Array2<f64>,
}
pub struct Phi {
    pub f: Array2<f64>,
    delta: Array2<f64>,
    pub params: Params,
}

pub struct Psi {
    pub f: Array2<f64>,
}


impl Field for Temperatura {
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}
impl Field for Phi{
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}
impl Field for Psi{
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}





unsafe fn conv_mul(v:&impl Field ,t:&impl Field,index:(usize,usize))->f64{
    return 2.0/std::f64::consts::PI*(v.dx(index)*t.dy(index)-v.dy(index)*t.dx(index))
}

impl Temperatura {
    pub fn new(nx: usize, ny: usize) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
            delta: Array2::<f64>::zeros((nx, ny)),
        };
    }
    pub fn step(&mut self, psi: &Psi,pr:f64, dt: f64) {
        unsafe {
            for i in 1..NX - 1 {
                for j in 1..NY - 1 {
                    let mut tmp = self.lap((i, j));
                    tmp += conv_mul(psi,self,(i, j));
                    tmp /= H * H;
                    *self.delta.uget_mut((i, j)) = tmp * dt;
                }
            }
        }
        self.f += &self.delta;
    }
}
impl Phi {
    pub fn new(nx: usize, ny: usize, params: Params) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
            delta: Array2::<f64>::zeros((nx, ny)),
            params: params,
        };
    }
    pub fn step(&mut self, psi: &Psi, t: &Temperatura,conc:&Concentration,  dt: f64,g:f64) {
        //let mut delta = Array2::<f64>::zeros((NX,NY));
        //let shape= delta.dim();
        unsafe {
            for i in 1..NX - 1 {
                for j in 1..NY - 1 {
                    let archim = self.params.rel*t.dx((i, j)) - self.params.params_c.bm*conc.dx((i, j)) ;
                    let mut tmp = self.params.pr*self.lap((i, j)) ;
                    tmp += 4.0/3.0*conv_mul(psi,self,(i, j));
                    tmp /= H * H;
                    tmp -= self.params.pr*self.f.uget((i,j))*std::f64::consts::PI.powi(2)/4.0;
                    tmp += self.params.pr*4.0/std::f64::consts::PI*g /H*archim;
                    //tmp -= self.params.pr * self.params.rel_c2 * (conc2.dx((i, j))) / H;
                    *self.delta.uget_mut((i, j)) = tmp * dt;
                }
            }
        }
        self.f += &self.delta;
    }


}
impl Psi {
    pub fn new(nx: usize, ny: usize) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
        };
    }
    pub fn step(&mut self, phi: &Phi) {
        poisson_relax(&phi.f, &mut self.f, H,PEREODIC);
    }
}

struct Q{
    pub f: Array2<f64>,
    
}
impl Q{
    pub fn new(nx: usize, ny: usize) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
        };
    }
}
impl Field for Q{
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}




pub struct Concentration {
    pub f: Array2<f64>,
    qew: Q,
    qsn: Q,
    vx: Array2<f64>,
    vy: Array2<f64>,
    pub params: ParamsConc,
}

impl Field for Concentration{
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}




impl Concentration {
    pub fn new(nx: usize, ny: usize, params: ParamsConc) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
            qew: Q::new(nx, ny),
            qsn: Q::new(nx, ny),
            vx: Array2::<f64>::zeros((nx, ny)),
            vy: Array2::<f64>::zeros((nx, ny)),
            params: params,
        };
    }
    fn calc_v(&mut self, psi: &Psi) {
        unsafe {
            for i in 1..NX {
                let mut k = 0;
                self.vx[[i, k]] = psi.f[[i, k + 1]] + psi.f[[i - 1, k + 1]];

                for k in 1..NY - 1 {
                    self.vx[[i, k]] = psi.dy((i, k)) + psi.dy((i - 1, k));
                }
                k = NY - 1;
                self.vx[[i, k]] = -psi.f[[i, k - 1]] - psi.f[[i - 1, k - 1]];
            }
            self.vx /= 2.0;

            for i in 1..NX - 1 {
                for k in 1..NY {
                    self.vy[[i, k]] = psi.dx((i, k)) + psi.dx((i, k - 1));
                }
            }

            for k in 1..NY {
                {
                    let i = 0;
                    self.vy[[i, k]] = psi.f[[i + 1, k]] + psi.f[[i + 1, k - 1]];
                }
                let i = NX - 1;
                self.vy[[i, k]] = -psi.f[[i - 1, k]] - psi.f[[i - 1, k - 1]];
            }
            self.vy /= 2.0;
           
        }
    }
    pub fn step(&mut self, psi: &Psi,temp:&Temperatura, dt: f64,g:f64) {
        let le = self.params.le;
        let sor= self.params.sor1;
        let l_sed = self.params.l_sed;
       //println!("{:?}",le/l);
        unsafe {
            self.calc_v(&psi);

            for i in 1..NX {
                for k in 0..NY {
                    let mut tmp = 2.0/std::f64::consts::PI*(self.vx[[i, k]] ) * self.mx_b((i, k));
                    tmp += -le * self.dx_b((i, k));
                    tmp += -le * sor*temp.dx_b((i, k))*self.mx_b((i,k));
                    self.qew.f[[i, k]] = tmp / H;
                }
            }
            for i in 0..NX {
                for k in 1..NY {
                    let mut tmp = (2.0/std::f64::consts::PI*(-self.vy[[i, k]]) -le/l_sed*H*g ) * self.my_b((i, k));
                    tmp += -le * (self.dy_b((i, k)));
                    tmp += -le * sor*temp.dy_b((i, k))*self.my_b((i,k)) ;
                    self.qsn.f[[i, k]] = tmp / H;
                }
            }
            /////////////////////////////////////////////////

            for i in 1..NX - 1 {
                for k in 1..NY - 1 {
                    let q = -(self.qew.dx_f((i, k)) + self.qsn.dy_f((i, k)));
                    self.f[[i, k]] += dt / H * q;
                }
            }
            for i in 1..NX - 1 {
                let k = 0;
                self.f[[i, k]] += dt / H * (-self.qew.dx_f((i, k)) - 2.0 * self.qsn.f[[i, k + 1]]);

                let k = NY - 1;
                self.f[[i, k]] += dt / H * (-self.qew.dx_f((i, k)) + 2.0 * self.qsn.f[[i, k]]);
            }
            if PEREODIC{
                for k in 0..NY {
                    self.f[[0, k]] = self.f[[NX - 2, k]];
                    self.f[[NX - 1, k]] = self.f[[1, k]];
                }
            }else{
                for k in 1..NY - 1 {
                    let i= 0;
                    self.f[[i, k]] += dt / H * (-self.qsn.dy_f((i, k)) - 2.0 * self.qew.f[[i+1, k]]);
    
                    let i = NX - 1;
                    self.f[[i, k]] += dt / H * (-self.qsn.dy_f((i, k)) +2.0 * self.qew.f[[i, k]]);
                }
                self.f[[0, 0]] += dt / H * 2.0 *( -self.qew.f[[1, 0]]- self.qsn.f[[0, 1]]);
                self.f[[0, NY-1]] +=  dt / H * 2.0 *( -self.qew.f[[1, NY-1]]+ self.qsn.f[[0, NY-1]]);
                self.f[[NX-1, 0]] += dt / H * 2.0 *( self.qew.f[[NX-1, 0]]- self.qsn.f[[NX-1, 1]]);
                self.f[[NX-1, NY-1]] += dt / H * 2.0 *( self.qew.f[[NX-1, NY-1]]+ self.qsn.f[[NX-1, NY-1]]);
            }
        }
    }
}


