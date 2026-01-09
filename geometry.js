//geometric library created by yrty2
//require geometric algebra in mathematics.js
//処理速度は考えられていない。(でもnewはそんなにコストにならないらしい)

//coordinate systems
class cartesian2D{
    constructor(x,y){
        this.x=x;
        this.y=y;
    }
    get length(){
        return Math.hypot(this.x,this.y);
    }
    get arg(){
        return Math.atan2(this.y,this.x);
    }
    add(v){
        return new cartesian2D(this.x+v.x,this.y+v.y);
    }
    sub(v){
        return new cartesian2D(this.x-v.x,this.y-v.y);
    }
    get poler(){
        return new spherical2D(this.length,this.arg);
    }
    scale(x){
        return new cartesian2D(this.x*x,this.y*x);
    }
    get vector(){
        const a=new Float64Array(2);
        a[0]=this.x;
        a[1]=this.y;
        return a;
    }
}
class cartesian{
    constructor(x,y,z){
        this.x=x;
        this.y=y;
        this.z=z;
    }
    scale(x){
        return new cartesian(this.x*x,this.y*x,this.z*x);
    }
    mul(cart){
        return new cartesian(this.x*cart.x,this.y*cart.y,this.z*cart.z);
    }
    normalize(){
        return this.normal;
    }
    get normal(){
        return this.scale(1/this.length);
    }
    get length(){
        return Math.sqrt(this.x*this.x+this.y*this.y+this.z*this.z);
    }
    get vector(){
        return [this.x,this.y,this.z];
    }
}
class cartesian4D{
    constructor(x,y,z,w){
        this.x=x;
        this.y=y;
        this.z=z;
        this.w=w;
    }
    scale(x){
        return new cartesian4D(this.x*x,
        this.y*x,
        this.z*x,
        this.w*x)
    }
    get length(){
        return Math.sqrt(this.x*this.x+this.y*this.y+this.z*this.z+this.w*this.w);
    }
}
class spherical2D{
    constructor(radius,theta){
        this.radius=radius;
        this.theta=theta;
    }
    convertCartesian(){
        return new cartesian2D(this.radius*Math.cos(this.theta),this.radius*Math.sin(this.theta));
    }
    get cartesian(){
        return this.convertCartesian();
    }
}
class polerSpherical{
    constructor(r,p,s){
        this.radius=r;
        this.theta=p;
        this.phi=s;
    }
    get cartesian(){
        //極はy軸とする。
        return (new cartesian(Math.cos(this.theta)*Math.sin(this.phi),Math.cos(this.phi),Math.sin(this.theta)*Math.sin(this.phi))).scale(this.radius);
    }
    translate(x,y){
        this.phi+=y;
        this.theta+=x;
    }
    get stereographic(){
        projection.stereographic(this);
    }
}
class spherical{
    //極が存在しない球面座標系
    constructor(r,v){
        this.r=r;
        this.x=v.x;
        this.y=v.y;
        this.z=v.z;
    }
    copy(){
        return new spherical(this.r,this.cartesian);
    }
    get clifford(){
        this.cliff=[0,this.x,this.y,this.z,0,0,0,0];
    }
    translate(x,y){
        //x=rで2piに相当
        const s=Math.hypot(x,y);
        if(s!=0){
        const rot=[Math.cos(Math.PI*s/this.r),0,0,0,0,y*Math.sin(Math.PI*s/this.r)/s,-x*Math.sin(Math.PI*s/this.r)/s,0];
        const v=clifford.rotate3D(new vector(this.x,this.y,this.z),rot);
        this.x=v.x;
        this.y=v.y;
        this.z=v.z;
        }
    }
    rotate(rot){
        const v=clifford.rotate3D(new vector(this.x,this.y,this.z),rot);
        this.x=v.x;
        this.y=v.y;
        this.z=v.z;
    }
    puremove(x,y){
        //x=rで2piに相当
        const s=Math.hypot(x,y);
        if(s!=0){
        const rot=[Math.cos(s/2),0,0,0,0,y*Math.sin(s/2)/s,-x*Math.sin(s/2)/s,0];
        const v=clifford.rotate3D(new vector(this.x,this.y,this.z),rot);
        this.x=v.x;
        this.y=v.y;
        this.z=v.z;
        }
    }
    rotor(p){
        //p,qの回転面(球面では単に軸)
        const c=clifford.product3D([0,this.x,this.y,this.z,0,0,0,0],[0,p.x,p.y,p.z,0,0,0,0]);
        const s=Math.sqrt(c[4]*c[4]+c[5]*c[5]+c[6]*c[6]);
        return [c[4]/s,c[5]/s,c[6]/s];
    }
    arg(p){
        //p,qの角度
        const l=Math.acos((this.x*p.x+this.y*p.y+this.z*p.z)/(this.r*p.r))
        if(isNaN(l)){
            return 0;
        }
        return l;
    }
    length(p){
        const l=Math.acos((this.x*p.x+this.y*p.y+this.z*p.z)/(this.r*p.r))*this.r;
        if(isNaN(l)){
            return 0;
        }
        return l;
    }
    translateBack(x,y){
        //x=rで2piに相当
        const s=Math.hypot(x,y);
        const rot=[Math.cos(Math.PI*s/this.r),0,0,0,0,y*Math.sin(Math.PI*s/this.r)/s,x*Math.sin(Math.PI*s/this.r)/s,0];
        const v=clifford.rotate3D(new vector(this.x,this.y,this.z),rot);
        return new spherical(this.r,v);
    }
    get cartesian(){
        return new cartesian(this.x,this.y,this.z);
    }
}
class spherical4D{
    constructor(r,v){
        this.r=r;
        this.x=v.x;
        this.y=v.y;
        this.z=v.z;
        this.w=v.w;
    }
    copy(){
        return new spherical4D(this.r,this.cartesian)
    }
    arg(p){
        //p,qの角度
        const l=Math.acos((this.x*p.x+this.y*p.y+this.z*p.z+this.w*p.w)/(this.r*p.r));
        if(isNaN(l)){
            return 0;
        }
        return l;
    }
    rotor(p){
        //p,qの回転面(球面では単に軸)
        const c=clifford.product4D([0,this.x,this.y,this.z,this.w,0,0,0,0,0,0,0,0,0,0,0],[0,p.x,p.y,p.z,p.w,0,0,0,0,0,0,0,0,0,0,0]);
        const s=Math.sqrt(c[5]*c[5]+c[6]*c[6]+c[7]*c[7]+c[8]*c[8]+c[9]*c[9]+c[10]*c[10]);
        if(s>0){
        return [c[5]/s,c[6]/s,c[7]/s,c[8]/s,c[9]/s,c[10]/s];
        }
        return -1;
    }
    get vector(){
        return [this.x,this.y,this.z,this.w];
    }
    cliffordrotate(rot){
        const v=clifford.rotate4D(this.vector,rot);
        this.x=v[0];
        this.y=v[1];
        this.z=v[2];
        this.w=v[3];
    }
    length(p){
        const l=this.arg(p)*this.r;
        if(isNaN(l)){
            return 0;
        }
        return l;
    }
    translate(x,y,z){
        const s=Math.sqrt(x*x+y*y+z*z);
        if(s>0){
        const theta=Math.PI*(s/this.r);
        const ss=Math.sin(theta)/s;
        //x=rで2piに相当
        const v=clifford.rotate4D([this.x,this.y,this.z,this.w],[Math.cos(theta),0,0,0,0,0,0,0,-x*ss,y*ss,-z*ss,0,0,0,0,0]);
        this.x=v[0];
        this.y=v[1];
        this.z=v[2];
        this.w=v[3];
        }
    }
    translateBack(x,y,z){
        //x=rで2piに相当
        const s=Math.sqrt(x*x+y*y+z*z);
        if(s>0){
        const theta=Math.PI*(s/this.r);
        const ss=Math.sin(theta)/s;
        //x=rで2piに相当
        const v=clifford.rotate4D([this.x,this.y,this.z,this.w],[Math.cos(theta),0,0,0,0,0,0,0,-x*ss,y*ss,-z*ss,0,0,0,0,0]);
        return new spherical4D(this.r,new cartesian4D(v[0],v[1],v[2],v[3]));
        }
        return this.copy();
    }
    rotate(xy,yz,zx){
        const s=Math.sqrt(xy*xy+yz*yz+zx*zx);
        if(s>0){
        const ss=Math.sin(s)/s;
        const v=clifford.rotate4D([this.x,this.y,this.z,this.w],[Math.cos(s),0,0,0,0,xy*ss,yz*ss,zx*ss,0,0,0,0,0,0,0,0]);
        this.x=v[0];
        this.y=v[1];
        this.z=v[2];
        this.w=v[3];
        }
    }
    get cartesian(){
        return new cartesian4D(this.x,this.y,this.z,this.w);
    }
}
class hyperbolic{
    //双曲面
    constructor(x,y){
        if(x instanceof Float32Array){
            this.z=x;
        }else{
        //use poincare disk
        const r=Math.hypot(x,y);
        if(r==0){
            this.z=c32.zero();
        }else{
        const hypt=Math.tanh(r/2)/r;
        this.z=c32.const(hypt*x,hypt*y);
        }
        }
    }
    translate(x,y){
        this.T(new hyperbolic(x,y));
    }
    translateBack(x,y){
        return this.TB(new hyperbolic(x,y));
    }
    T(a){
        this.z=this.TB(a.z);
    }
    TB(a){
        return c32.quot(c32.sum(a,this.z),c32.sum(c32.mul(this.z,c32.conjugate(a)),c32.const(1,0)));
    }
    R(theta){
        this.z=this.RB(theta);
    }
    RB(theta){
        return c32.mul(this.z,c32.poler(1,theta));
    }
    Ra(a,theta){
        this.z=Rab(a,theta);
    }
    RaB(a,theta){
        return this.TB(c32.mul(this.TB(this.z,c32.neg(a.z)),c32.poler(1,theta)),a.z);
    }
    get length(){
        return 2*Math.atanh(c32.abs(this.z)[0]);
    }
    dist(a){
        return 2*Math.atanh(c32.abs(this.TB(c32.neg(a.z)))[0]);
    }
    get arg(){
        return c32.arg(this.z)[0];
    }
    angle(a){
        return c32.arg(hyperbolicGeometry.translate(a.z,c32.neg(this.z)))[0];
    }
    geodesic(a,d){
        if(!d){
            d=12;
        }
        //aまでの直線
        const points=[];
        for(let k=0; k<=d; ++k){
            const theta=this.angle(a);
            const t=this.dist(a)*k/d;
            points.push(hyperbolicGeometry.translate(c32.poler(Math.tanh(t/2),theta),this.z));
        }
        return points;
    }
    plot(canvas,context){
        context.beginPath();
        context.arc((canvas.height*this.x+canvas.width)/2,canvas.height*(-this.y+1)/2,5,0,2*Math.PI);
        context.fill();
        context.closePath();
    }
    get vector(){
        return this.z;
    }
    get klein(){
        //ベルトラミ・クラインモデル
        const abs=c32.abs(this.z)[0];
        return c32.prod(c32.prod(this.z,2),1/(1+abs*abs));
    }
    get upperhalf(){
        //上半平面
        const z=c32.const(this.y,-this.x);   
        return c32.sub(c32.mul(c32.quot(c32.sum(c32.one(),z),c32.sub(c32.one(),z)),c32.imag(1)),c32.imag(1));
    }
    disk(f){
        //0でクライン、1でポアンカレ
        const abs=c32.abs(this.z)[0];
        return c32.prod(c32.prod(this.z,2),1/(1+f+abs*abs*(1-f)));
    }
    scale(a){
        const abs=c32.abs(this.z)[0];
        this.z=c32.prod(this.z,Math.tanh(a*Math.atanh(abs))/abs);
    }
    SB(a){
        const abs=c32.abs(this.z)[0];
        return c32.prod(this.z,Math.tanh(a*Math.atanh(abs))/abs);
    }
    get x(){
        return this.z[0];
    }
    get y(){
        return this.z[1];
    }
}
const hyperbolicGeometry={
    distance(p,q){
        return 2*Math.atanh(c32.abs(this.translate(p,c32.neg(q)))[0]);
    },
    arg(p){
        return c32.arg(p)[0];
    },
    midpoint(p,q){
        //pもqも複素数であるとする。
        return this.translate(c32.poler(Math.tanh(this.distance(p,q)/4),this.arg(this.translate(p,c32.neg(q)))),q);
    },
    translate(z,a){
        return c32.quot(c32.sum(a,z),c32.sum(c32.mul(z,c32.conjugate(a)),c32.one()));
    },
    //鏡映変換
    reflection(p,q,point){
        const m=this.midpoint(p,q);
        return this.translate(c32.neg(point),this.scale(m,2));
    },
    reflectionPolygon(p,q,H){
        const m=this.midpoint(p,q);
        //中点m,多角形H
        const res=[];
        for(const h of H){
            if(h instanceof hyperbolic){
                res.push(new hyperbolic(this.translate(c32.neg(h.z),this.scale(m,2))));
            }else{
                res.push(this.translate(c32.neg(h),this.scale(m,2)));
            }
        }
        return res;
    },
    scale(p,a){
        const abs=c32.abs(p)[0];
        return c32.prod(p,Math.tanh(a*Math.atanh(abs))/abs);
    }
}
class topology{
}
const projection={
    //射影
    orthogonal(cart){
        if(cart.z>0){
        return new cartesian2D(cart.x,cart.y);
        }
    },
    perspective(cart){
        if(cart.z>0){
            return new cartesian2D(cart.x/cart.z,cart.y/cart.z);
        }
    },
    stereographic(pole){
        if(pole.constructor==polerSpherical){
        return (new cartesian2D(Math.cos(pole.theta)*Math.sin(pole.phi),
                Math.sin(pole.theta)*Math.sin(pole.phi))).scale(pole.radius/(1-Math.cos(pole.phi)));
        }
        return (new cartesian2D(pole.x,pole.y)).scale(1/(1-pole.z/pole.r));
    },
    stereographic3D(pole){
        return (new cartesian(pole.x,pole.y,pole.z)).scale(1/(1-pole.w/pole.r));
    },
    poincareDisk(hyp){
        return this.orthogonal(hyp);
        //return (new cartesian2D(hyp.x,-hyp.y)).scale(1/(1+Math.sqrt(1+hyp.x*hyp.x+hyp.y*hyp.y)));
    }
}
function Cl(hyperbolic,imaginary){
    let V=Array(hyperbolic).fill(1);
    V.push(...Array(imaginary).fill(-1));
    //v is like [1,1,-1] [-1,-1,-1]
    const n=hyperbolic+imaginary;
    let cl=Array(n);
    for(let k=0; k<n; k++){
        cl[k]=k;
    }
    cl=maths.power(cl);
    //返すのは数式
    function Clmul(u,v){
        //u,v->[0,1,2] これは基底の積。最終的に昇順に。
        let a=[...u,...v];//[0,1,1,2],[1]^2は？
        let h=1;
        //入れ替えソート(符号反転が行われる。)
        while(true){
            let zyun=true;
            let hold=0;
            for(let k=0; k<a.length; ++k){
                if(hold<=a[k]){
                }else{
                    zyun=false;
                    break;
                }
                hold=a[k];
            }
            if(zyun){
                break;
            }
            //ここに処理
            for(let k=1; k<a.length; ++k){
                if(a[k-1]>a[k]){
                    const holder=a[k-1];
                    a[k-1]=a[k];
                    a[k]=holder;
                    h*=(-1);
                }
            }
        }
        //2乗項を探す。(場合によっては符号反転が行われる)
        for(let k=1; k<a.length; ++k){
            if(a[k-1]==a[k]){
                h*=V[a[k]];
                a=[...a.slice(0,k-1),...a.slice(k+1,a.ength)]
                k--;
            }
        }
        return [a,h];
    }
    const tapes=Array(cl.length).fill("");
    //冪集合の積
    let tape="return [";
    for(let i=0; i<cl.length; ++i){
    for(let j=0; j<cl.length; ++j){
        const a=Clmul(cl[i],cl[j]);
        const id=cl.findIndex(e=>e.join()==a[0].join());
        if(id!=-1){
            let hugou="+";
            if(a[1]==-1){
                hugou="-";
            }else if(tapes[id].length==0){
                hugou="";
            }
            tapes[id]+=`${hugou}p[${i}]*q[${j}]`;
        }else{
            console.warn("おい！おかしいぞ！");
        }
    }
    }
    for(let k=0; k<tapes.length; ++k){
        tape+=tapes[k];
        if(k+1<tapes.length){
            tape+=",";
        }
    }
    return tape+"]";
}
const c32={
    const(a,b){
        const c=new Float32Array(2);
        c[0]=a;
        c[1]=b;
        return c;
    },
    one(){
        const c=new Float32Array(2);
        c[0]=1;
        return c;
    },
    real(a){
        const c=new Float32Array(2);
        c[0]=a;
        return c;
    },
    imag(a){
        const c=new Float32Array(2);
        c[1]=a;
        return c;
    },
    zero(){
        return new Float32Array(2);
    },
    neg(z){
        const c=new Float32Array(2);
        c[0]=-z[0];
        c[1]=-z[1];
        return c;
    },
    poler(radius,theta){
        const c=new Float32Array(2);
        c[0]=radius*Math.cos(theta);
        c[1]=radius*Math.sin(theta);
        return c;
    },
    prod(z,x){
        return this.const(z[0]*x,z[1]*x);
    },
    exp(z){
        const r=Math.exp(z[0]);
        return this.const(r*Math.cos(z[1]),r*Math.sin(z[1]));
    },
    mul(z,w){
        if(!Array.isArray(z) && !Array.isArray(w)){
            const c=new Float32Array(2);
            c[0]=z[0]*w[0]-z[1]*w[1];
            c[1]=z[0]*w[1]+z[1]*w[0];
            return c;
        }
        if(!Array.isArray(z) && Array.isArray(w)){
            return c32v.mul(w,z);
        }
        if(Array.isArray(z) && !Array.isArray(w)){
            return c32v.mul(z,w);
        }
        return c32v.cross(z,w);
    },
    sum(z,w){
        if(Array.isArray(z)){
            return c32v.sum(z,w);
        }
        const c=new Float32Array(2);
        c[0]=z[0]+w[0];
        c[1]=z[1]+w[1];
        return c;
    },
    sub(z,w){
        if(Array.isArray(z)){
            return c32v.sub(z,w);
        }
        const c=new Float32Array(2);
        c[0]=z[0]-w[0];
        c[1]=z[1]-w[1];
        return c;
    },
    abs(z){
        if(Array.isArray(z)){
            return c32v.length(z);
        }
        const c=new Float32Array(2);
        c[0]=Math.sqrt(z[0]*z[0]+z[1]*z[1]);
        return c;
    },
    normalize(z){
        return this.prod(z,1/Math.sqrt(z[0]*z[0]+z[1]*z[1]));
    },
    conjugate(z){
        const c=new Float32Array(2);
        c[0]=z[0];
        c[1]=-z[1];
        return c;
    },
    quot(z,w){
        if(w[1]==0){
            const c=new Float32Array(2);
            c[0]=z[0]/w[0];
            c[1]=z[1]/w[0];
            return c;
        }
        return this.prod(this.mul(z,this.conjugate(w)),1/(w[0]*w[0]+w[1]*w[1]));
    },
    arg(z){
        return this.const(Math.atan2(z[1],z[0]),0);
    },
    log(z){
        return this.const(Math.log(z[0]*z[0]+z[1]*z[1])/2,Math.atan2(z[1],z[0]));
    },
    pow(z,w){
        if(w[1]==0){
            const c=new Float32Array(2);
            const theta=w[0]*Math.atan2(z[1],z[0]);
            const r=Math.pow(z[0]*z[0]+z[1]*z[1],w[0]/2);
            c[0]=r*Math.cos(theta);
            c[1]=r*Math.sin(theta);
            return c;
        }else{
            //複雑な式
            const theta=Math.atan2(z[1],z[0]);
            const lnr=Math.log(z[0]*z[0]+z[1]*z[1])/2;
            const r=Math.exp(w[0]*lnr-w[1]*theta);
            const phi=w[0]*theta+w[1]*lnr;
            return this.const(r*Math.cos(phi),r*Math.sin(phi));
        }
    },
    fact(z){
        if(z[0]>=0 && z[0]%1==0 && z[1]==0){
            //自然数階乗
            let res=1;
            for(let k=1; k<=z[0]; ++k){
                res*=k;
            }
            return this.const(res,0);
        }else{
        if(z[0]<-1){
            const zaddone=this.const(z[0]+1,z[1]);
            return this.quot(this.fact(zaddone),zaddone);
        }else{
        let res=c32.const(0,0);
        for(let k=1; k<=10000; ++k){
            const rz=Math.exp(-k/100);
            const fz=this.pow(c32.const(k/100,0),z);
            res[0]+=rz*fz[0];
            res[1]+=rz*fz[1];
        }
            res[0]*=1/100;
            res[1]*=1/100;
        return res;
        }
        }
    },
    factn(n,x){
        if(x[1]!=0 || x[0]%1!=0){
            console.error("多重階乗は複素数または実数に拡張されていません！");
            return;
        }
        let res=1;
        if(x[0]<0){
            if((-x[0])%n==0){
                /*(x+2)!!/(x+2)=x!! (-N*n)!^n=NaN*/
                return c32.const(NaN,NaN);
            }else{
                return c32.prod(c32.factn(n,c32.const(x[0]+n,0)),1/(x[0]+n));
            }
        }else{
        for(let k=x[0]; k>0; k-=n){
            res*=k;
        }
        return c32.const(res,0);
        }
    },
    print(z){
        let operator="+";
        let x=Math.round(z[0]*1000)/1000;
        let y=Math.round(z[1]*1000)/1000;
        if(x==0){
            x="";
            if(operator=="+"){
            operator="";
            }
        }
        if(y<0){
            operator="";
        }
        if(y==0){
            return `${x}`;
        }
        if(y==1){
            y="";
        }
        if(y==-1){
            y="-";
        }
        return `${x}${operator}${y}i`;
    },
    re(z){
        return c32.const(z[0],0);
    },
    im(z){
        return c32.const(0,z[1]);
    },
    floor(z){
        return c32.const(Math.floor(z[0]),Math.floor(z[1]));
    },
    mod(z,w){
        if(z[1]==0 && w[1]==0){
            return c32.const(z[0]-w[0]*Math.floor(z[0]/w[0]),0);
        }
        return c32.sub(z,c32.mul(w,c32.floor(c32.quot(z,w))));
    },
    sin(z){
        if(z[1]==0){
            return c32.const(Math.sin(z[0]),0);
        }
        const ez=c32.exp(c32.const(-z[1],z[0]));
        const mez=c32.exp(c32.const(z[1],-z[0]));
        return c32.quot(c32.sub(ez,mez),c32.const(0,2));
    },
    cos(z){
        if(z[1]==0){
            return c32.const(Math.cos(z[0]),0);
        }
        const ez=c32.exp(c32.const(-z[1],z[0]));
        const mez=c32.exp(c32.const(z[1],-z[0]));
        return c32.prod(c32.sum(ez,mez),1/2);
    },
    tan(z){
        const ez=c32.exp(c32.const(-2*z[1],2*z[0]));
        const ezu=c32.const(ez[0]-1,ez[1]);
        const ezd=c32.const(ez[0]+1,ez[1]);
        return c32.mul(c32.const(0,-1),c32.quot(ezu,ezd));
    },
    sqrt(z){
        if(z[1]==0 && z[0]>=0){
            return c32.const(Math.sqrt(z[0]),0);
        }
        const r=Math.pow(z[0]*z[0]+z[1]*z[1],0.25);
        const theta=Math.atan2(z[1],z[0])/2;
        return c32.const(r*Math.cos(theta),r*Math.sin(theta));
    },
    tetration(z,w){
        //wは整数である必要がある。
        if(w[0]==-1){
            return new Float32Array(2);
        }
        let res=c32.const(1,0);
        for(let k=1; k<=w[0]; ++k){
            res=c32.pow(z,res);
        }
        return res;
    },
    int(f,C){
        //不定積分しないよ。
        //シンプソン法の方が早い...？不明なので区分求積法を使う。
        if(C[0][1]==0 && C[1][1]==0){
        const I=100000;
        let res=0;
        for(let k=1-C[0][0]*I; k<=C[1][0]*I; ++k){
            res+=f(k/I);
        }
        return new Float32Array([res/I,0]);
        }
        //経路積分は経路が与えられないといけない。
        console.warn("経路積分は未実装");
    },
    der(f){
        //fは関数である必要がある。
        //浮遊小数点の限界を感じる。
        const h=0.01;
        return z=>c32.prod(c32.sub(f(c32.const(z[0]+h,z[1])),f(z)),1/h);
    },
    plot(z,canvas,context){
        context.beginPath();
        context.arc((canvas.height*z[0]+canvas.width)/2,canvas.height*(-z[1]+1)/2,5,0,2*Math.PI);
        context.fill();
        context.closePath();
    },
    line(pointList,canvas,context){
        context.beginPath();
        for(const p of pointList){
        context.lineTo((canvas.height*p[0]+canvas.width)/2,canvas.height*(-p[1]+1)/2);
        }
        context.stroke();
        context.closePath();
    }
}