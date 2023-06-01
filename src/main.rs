const NY:usize=33;
const L:f64=2.0;
const NX:usize=(NY-1)*(L as usize)+1;

const H:f64=L/((NX-1) as f64);
const PEREODIC:bool = true;

use ndarray_stats::QuantileExt;
pub mod equations;
use equations::{Temperatura,Psi,Phi,Params,Concentration};
use csv::Writer;
use serde::{Serialize};
use std::fs;


use pdifflib::io::{write_stage};
use pdifflib::system::System;
use toml;


struct S{
    temp:Temperatura,
    psi:Psi,
    phi:Phi,
    conc:Concentration,
    params:Params
}



impl System for S{
    fn next_step(&mut self,_dt:f64,time:f64){
        self.phi.step(&self.psi,&self.temp,&self.conc,_dt,time);
        self.psi.step(&self.phi);
        self.temp.step(&self.psi,self.params.pr,_dt);
        self.conc.step(&self.psi,&self.temp,_dt);
        self.boundary_condition();
    }
    fn write_mat(&mut self, time:f64){
        let stage  = vec![&self.phi.f,&self.psi.f,&self.temp.f,&self.conc.f];
        let header = vec![String::from("x"),String::from("y"),String::from("phi"),String::from("psi"),String::from("T"),String::from("C")];

        write_stage(stage,header, String::from(format!("res/stage_t={:05}",time)),H)

    }
    fn boundary_condition(&mut self){
        if PEREODIC{

            for i in 0..NY{

                self.phi.f[[0,i]]=self.phi.f[[NX-2,i]];
                self.phi.f[[NX-1,i]]=self.phi.f[[1,i]];

                self.temp.f[[0,i]]=self.temp.f[[NX-2,i]];
                self.temp.f[[NX-1,i]]=self.temp.f[[1,i]];

            }
        }else{
            for i in 1..NY-1{
                self.phi.f[[0,i]]=-self.psi.f[[1,i]]/(H*H)*2.0;
                self.phi.f[[NX-1,i]]=-self.psi.f[[NX-2,i]]/(H*H)*2.0;
            }
        }

        for i in 0..NX{
            self.phi.f[[i,0]]=-self.psi.f[[i,1]]/(H*H)*2.0;
            self.phi.f[[i,NY-1]]=-self.psi.f[[i,NY-2]]/(H*H)*2.0;
        }

    }
    fn log_params(&self,wrt:&mut Writer<fs::File>,time:f64){
        let mut nu  =0.0;
        let shape = self.temp.f.dim();
        for i in 1..shape.0-1{
            nu+= self.temp.f[[i,shape.1-1]]-self.temp.f[[i,shape.1-2]]
        }
        nu/=H*(shape.0 as f64)-2.0;

       let row = Row{
            t:time,
            psi_m:*(self.psi.f.max().unwrap()),
            psi_l:(self.psi.f[[NX/4,NY/2]]),
            nu:nu,
       };
        wrt.serialize(&row).unwrap();

    }
    fn get_max_time(&self)->f64{
        return self.params.time;
    }
    fn get_DT(&self)->f64{
        return H*H/5.0/10.;
    }
    fn initial_condition(&mut self){
        let gama:f64= self.params.params_c.sor;
 
      for i in 0..NX{
        for j in 0..NY{

            self.temp.f[[i,j]] = 1.0-(j as f64)*H;
            let z = (j as f64)*H;
            let x= (i as f64)*H;
            let pi =std::f64::consts::PI;
            self.phi.f[[i,j]] =2.00;

            self.conc.f[[i,j]] =gama*(1.0-z) //gama*f64::exp(-z*gama)/(1.0-f64::exp(-gama));
        }
    }


    }
}

#[derive(Serialize)]
pub struct Row{
    t:f64,
    psi_m:f64,
    psi_l:f64,
    nu:f64,
}
fn read_params()->Params{
    let contents = fs::read_to_string("config.toml").unwrap() ;
    let params:Params = toml::from_str(&contents).unwrap();
    params
}

fn main() {
    println!("Hello, world!");
    //let mut dela =Array2::<f64>::zeros((NX,NY));
    let params = read_params();
    println!("{:?}",params);

    let mut system = S{
        temp:Temperatura::new(NX,NY),
        phi:Phi::new(NX,NY,params),
        psi:Psi::new(NX,NY),
        params:params,
        conc:Concentration::new(NX,NY,params.params_c)
    };
    system.solve()

}
