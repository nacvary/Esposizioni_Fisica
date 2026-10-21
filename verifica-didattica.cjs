// Independent checks for work integrals, free motion and charge assembly.
// Run: node verifica-didattica.cjs. Node is not needed to use the HTML.
const fs=require('fs'),vm=require('vm'),assert=require('assert'),path=require('path');
const html=fs.readFileSync(path.join(__dirname,'index.html'),'utf8');
const scripts=[...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m=>m[1]).join('\n');
const names=['fieldAt','workPathPoint','integrateWorkPath','motionAt','motionWindow','pairEnergy','assemblyData','ringPotential'];
let code='const EPS0=8.8542e-12,K=1/(4*Math.PI*EPS0);\n';
for(const name of names){const line=scripts.split('\n').find(l=>l.startsWith('function '+name+'('));assert(line,'Missing '+name);code+=line+'\n';}
code+='this.api={K,'+names.join(',')+'};';const sandbox={};vm.runInNewContext(code,sandbox);const A=sandbox.api;
let checks=0,maxIntegralRelativeError=0,maxEnergyRelativeDrift=0;
function ok(condition,message){checks++;assert(condition,message);}
function near(a,b,rtol=1e-9,atol=1e-11){ok(Number.isFinite(a)&&Number.isFinite(b)&&Math.abs(a-b)<=atol+rtol*Math.max(Math.abs(a),Math.abs(b)),`${a} != ${b}`);}
const dipole=[{q:1000,x:-100,y:0,z:0},{q:-1000,x:100,y:0,z:0}],single=[{q:1,x:0,y:0,z:0}];
const start=[-60,80],end=[60,80];
let seed=35713;function rand(){seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;}
// Use points and curves entirely above all sources, then independent endpoint potentials.
for(let j=0;j<60;j++){
  const charges=Array.from({length:3},()=>({q:rand()*2000-1000,x:rand()*240-120,y:-60-rand()*80,z:rand()*100-50}));
  const a=[-150+rand()*40,50+rand()*30],b=[110+rand()*40,50+rand()*30],bend=rand()*70;
  const exact=A.fieldAt(charges,...a,0).V-A.fieldAt(charges,...b,0).V;
  for(const curve of [0,bend]){const integral=A.integrateWorkPath(charges,a,b,curve);ok(integral.valid,'Safe path invalid');near(integral.integral,exact,1e-8,1e-4);maxIntegralRelativeError=Math.max(maxIntegralRelativeError,Math.abs(integral.integral-exact)/Math.max(1,Math.abs(exact)));}
}
const gamma=A.integrateWorkPath(dipole,start,end),delta=A.integrateWorkPath(dipole,start,end,70);
near(gamma.integral,delta.integral,1e-9,1e-3);near(A.integrateWorkPath(dipole,end,start).integral,-gamma.integral,1e-10,1e-4);
// Charge sign and units: pC -> C, µm -> m; independent radial antiderivative.
const radial=A.integrateWorkPath(single,[50,0],[150,0]);near(radial.integral,A.K*1e-9*(1/50e-6-1/150e-6),1e-10,1e-5);
near(1e-12*radial.integral,-(-1e-12*radial.integral));near(0*radial.integral,0);
const circle=A.integrateWorkPath(single,[100,0],[0,100],0,true),chord=A.integrateWorkPath(single,[100,0],[0,100]);ok(circle.valid&&chord.valid);near(circle.integral,0,0,1e-8);near(chord.integral,0,0,1e-8);
for(let t=0;t<=1;t+=.05){const {p,d}=A.workPathPoint([100,0],[0,100],0,t,true),f=A.fieldAt(single,...p);near(Math.hypot(...p),100,1e-12);near((f.Ex*d[0]+f.Ey*d[1])*1e-6,0,0,1e-8);}
ok(!A.integrateWorkPath(single,[-100,0],[100,0]).valid,'Crossing singularity must be rejected');
ok(!A.integrateWorkPath(single,[0,0],[100,0]).valid,'Singular endpoint must be rejected');
ok(A.integrateWorkPath([{q:0,x:0,y:0,z:0}],[-100,0],[100,0]).valid,'Zero source does not create singularity');
// Energy conservation and the work-energy theorem, including both charge signs and turning.
const base={q:1,mass:1,v0:0,voltage:100,distance:10,start:20};
for(let j=0;j<100;j++){
  const p={q:(rand()-.5)*20,mass:.1+rand()*99.9,v0:rand()*4-2,voltage:rand()*1000,distance:1+rand()*49,start:10+rand()*80};
  const win=A.motionWindow(p),initial=A.motionAt(p,0);ok(win.reachesEdge&&win.duration>0);
  for(let i=0;i<=30;i++){
    const f=A.motionAt(p,win.duration*i/30),expectedTotal=.5*p.mass*1e-6*p.v0**2+p.q*1e-9*p.voltage*(1-p.start/100);
    near(f.total,expectedTotal,2e-10,1e-17);near(f.kinetic-initial.kinetic,f.work,2e-10,1e-17);
    ok(f.x>=-1e-12&&f.x<=f.d+1e-12,'Trajectory exits before declared duration');
    maxEnergyRelativeDrift=Math.max(maxEnergyRelativeDrift,Math.abs(f.total-expectedTotal)/Math.max(1e-20,Math.abs(initial.U)+initial.kinetic));
  }
  const final=A.motionAt(p,win.duration);near(Math.min(Math.abs(final.x),Math.abs(final.x-final.d)),0,0,1e-10);
}
const stop=A.motionWindow({...base,q:0,v0:0});ok(!stop.reachesEdge);near(A.motionAt({...base,q:0},.3).x,.02);
const inertial={...base,voltage:0,v0:.3},iw=A.motionWindow(inertial);near(iw.duration,(.1-.02)/.3);near(A.motionAt(inertial,iw.duration).v,.3);
const turn={...base,start:50,v0:-.12},turnT=.12/A.motionAt(turn,0).a;ok(turnT<A.motionWindow(turn).duration);near(A.motionAt(turn,turnT).v,0);ok(A.motionAt(turn,turnT).x>0);ok(A.motionAt(turn,turnT*2).v>0);
near(A.motionAt({...base,mass:2},0).a,A.motionAt(base,0).a/2);near(A.motionAt({...base,q:-1},0).a,-A.motionAt(base,0).a);
// Assembly excludes self-energy and counts each pair once. Test true 3D distances.
for(let n=1;n<=6;n++){
  const charges=Array.from({length:n},(_,i)=>({q:(i%2?-1:1)*(i+1)*100,x:i*20-60,y:i*i*3,z:i*17}));let reference=0;
  for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){const x=charges[i],y=charges[j],distance=Math.hypot(x.x-y.x,x.y-y.y,x.z-y.z)*1e-6;reference+=A.K*x.q*y.q*1e-18/distance;}
  for(const reverse of [false,true]){let incremental=0;for(let step=0;step<=n;step++){const data=A.assemblyData(charges,reverse,step);ok(data.pairs.length===step*(step-1)/2);if(step)incremental+=data.increment;near(data.total,incremental,1e-12);if(step===n)near(data.total,reference,1e-12);} }
}
near(A.assemblyData(dipole,false,1).total,0);ok(A.assemblyData(dipole,false,2).total<0);ok(A.assemblyData([{q:1,x:0,y:0,z:0},{q:1,x:1,y:0,z:0}],false,2).total>0);
ok(!A.assemblyData([{q:1,x:0,y:0,z:0},{q:-1,x:0,y:0,z:0}],false,2).valid);near(A.pairEnergy({q:0,x:0,y:0,z:0},{q:1,x:0,y:0,z:0}),0);
for(const n of [3,12,60,120])for(const z of [-200,-60,0,60,200]){const v=A.ringPotential(n,z);near(v.discrete,v.continuous,1e-12,1e-8);}
const report={status:'PASS',checks,randomWorkConfigurations:60,randomMotionConfigurations:100,maxIntegralRelativeError,maxEnergyRelativeDrift,defaultWorkJ:gamma.integral*1e-12,closedLoopResidualJ:(gamma.integral-delta.integral)*1e-12,equipotentialArcIntegralV:circle.integral,assemblyChargeCounts:[1,2,3,4,5,6]};
fs.writeFileSync(path.join(__dirname,'risultati-didattica.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report,null,2));
