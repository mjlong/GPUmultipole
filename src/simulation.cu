#include "simulation.h" 
#define CHOP 0.7
#define NU2 0.55  //0.5<-->2.5; 0.55<-->2.45
extern __constant__ float wdspp[];

__global__ void fixsrc_sample(MemStruct pInfo, float width, int shift){
  int idr = blockDim.x * blockIdx.x + threadIdx.x;
  int id  = idr + shift;
  curandState state = pInfo.nInfo.rndState[idr];
  pInfo.nInfo.pos_x[id] =width*curand_uniform_double(&state);// 
  pInfo.nInfo.pos_y[id] =width*curand_uniform_double(&state);// 
  pInfo.nInfo.pos_z[id] =width*curand_uniform_double(&state);// 
  pInfo.nInfo.rndState[idr] = state;
}

__global__ void initialize(MemStruct pInfo,float width, int banksize,int shift, int seed){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x + shift;
  /* Each thread gets same seed, a different sequence number, no offset */
  //curand_init(id*7546861334684321478, id, id+14412078966483154, &(pInfo.nInfo.rndState[id]));
  curand_init(9798+seed, id, 0, &(pInfo.nInfo.rndState[id-shift]));
  neutron_sample(pInfo.nInfo, id, id-shift, width);
#if defined(__CTALLY)
  pInfo.tally.cnt[id-shift] = 0;
#if defined(__CTALLY2)
  pInfo.tally.cnt2[id-shift] = 0;
#endif
#endif 
#if defined(__MTALLY)
#if defined(__1D)
  pInfo.nInfo.imat[id] = int(floorf(pInfo.nInfo.pos_x[id]/wdspp[1]));
#endif
#if defined(__3D)
  pInfo.nInfo.imat[id] = ((int)floorf(pInfo.nInfo.pos_x[id]/wdspp[1]) + 
			  (int)floorf(pInfo.nInfo.pos_y[id]/wdspp[1])*(int)wdspp[5] + 
			  (int)floorf(pInfo.nInfo.pos_z[id]/wdspp[1])*(int)(wdspp[5]*wdspp[5]));
#endif
#endif
  pInfo.nInfo.live[id] = 1*(id<banksize);
}

__device__ void neutron_sample(NeutronInfoStruct nInfo, unsigned id,unsigned idr, float width){
  nInfo.live[id] = 1u;
  curandState state = nInfo.rndState[idr];
  //TODO: source sampling should take settings dependent on geometry
#if defined(__1D)
#if defined(__MTALLY)
  nInfo.pos_x[id] =width*curand_uniform_double(&state);
#else
  nInfo.pos_x[id] =width/(CHOP*PI)*asin(sin(PI*0.5*CHOP)*(1-2*curand_uniform_double(&state)))+width*0.5;//width*curand_uniform_double(&state);
#endif
#endif
#if defined(__3D)
  //nInfo.pos_x[id] =width/PI*acos(1-2*curand_uniform_double(&state));
  //width*curand_uniform_double(&state);// 
  //nInfo.pos_y[id] =width/PI*acos(1-2*curand_uniform_double(&state));
  //width*curand_uniform_double(&state);// 
  //nInfo.pos_z[id] =width/PI*acos(1-2*curand_uniform_double(&state));
  //width*curand_uniform_double(&state);// 
  //nInfo.pos_x[id] =width*curand_uniform_double(&state);
  //width/(CHOP*PI)*asin(sin(PI*0.5*CHOP)*\(1-2*curand_uniform_double(&state)))+width*0.5; 
  //nInfo.pos_y[id] =width*curand_uniform_double(&state);
  //width/(CHOP*PI)*asin(sin(PI*0.5*CHOP)*(1-2*curand_uniform_double(&state)))+width*0.5; 
  //nInfo.pos_z[id] =width*curand_uniform_double(&state);
  //width/(CHOP*PI)*asin(sin(PI*0.5*CHOP)*(1-2*curand_uniform_double(&state)))+width*0.5; 
  nInfo.pos_x[id] =width*0.5;
  nInfo.pos_y[id] =width*0.5;
  nInfo.pos_z[id] =width*0.5;
#endif
  nInfo.rndState[idr] = state;
#if defined(__WASTE)
  nInfo.energy[id] = STARTENE;
#endif
#if defined(__SCATTERPLOT)
  nInfo.energy[id] = nInfo.pos_z[id];
#endif
}

#if defined(__3D)
__device__ float intersectbox(float x, float y, float z, float a, float b, float c, float vx, float vy, float vz){
  //float t1,t2;
  /*t1 = max(min(-x/vx,(a-x)/vx),max(min(-y/vy,(b-y)/vy),min(-z/vz,(c-z)/vz)));*/
  return min(max(-x/vx,(a-x)/vx),min(max(-y/vy,(b-y)/vy),max(-z/vz,(c-z)/vz)));
  /*printf("\n id=%2d,from (%.10e,%.10e,%.10e) along (%.10e,%.10e,%.10e), to\n (%.10e,%.10e,%.10e) or (%.10e,%.10e,%.10e),\n",
	 blockDim.x * blockIdx.x + threadIdx.x,x,y,z,vx,vy,vz,
	 x+vx*t1,y+vy*t1,z+vz*t1, x+vx*t2,y+vy*t2,z+vz*t2);*/
  //return t2; //in this box case, it's certain t1<0<t2, and t2 is what's used
}


__device__ float product(float* n1, float* n2){
  return n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
}

__device__ void add(float *v1, float* v2, float multi){
  v1[0]+=v2[0]*multi;
  v1[1]+=v2[1]*multi;
  v1[2]+=v2[2]*multi;
}

__global__ void preview_live(MemStruct DeviceMem, int shift){
  int id = blockDim.x * blockIdx.x + threadIdx.x + shift;
  int live = DeviceMem.nInfo.live[id];
  if(34217==id) printf("  checking... live[34217]=%d\n",live);
  if(live>3) printf("  checking... id=%d,live=%d>3\n",id,live);
}

__device__ int calind(float x, float dx){
  int i=0;
  while(dx*i<x)
    i++;
  return i-1;
}

__global__ void history(MemStruct DeviceMem, unsigned num_src,int shift,unsigned banksize){
  //float a = wdspp[0];
  //float dx = wdspp[1];
  //float mfp = wdspp[2];
  //float Ps = 1-(wdspp[3]+wdspp[4]);
  //float Pc = Ps+wdspp[4]; //=1-wdspp[3]

  float n[3] = {0.0, 0.0,0.0};
  float v[3];

  float l;//sampled path
  float t;//length to boundary
  float s;

  int id = blockDim.x * blockIdx.x + threadIdx.x + shift;
  curandState localState = DeviceMem.nInfo.rndState[id-shift];
#if defined(__CTALLY2)||(__FTALLY2)||(__MTALLY)
  int nid = id+num_src;
#else
  int nid = int(curand_uniform_double(&localState)*banksize)+num_src;
#endif
  //if(100>id) {printf("  id=%d,nid=%d,rnd=%.5f\n",id,nid,curand_uniform_double(&localState));}
  //extern __shared__ unsigned blockTerminated[];

  CMPTYPE rnd;
  float x = DeviceMem.nInfo.pos_x[nid];
  float y = DeviceMem.nInfo.pos_y[nid];
  float z = DeviceMem.nInfo.pos_z[nid];
  float mu,phi; 

#if defined(__FTALLY)
  DeviceMem.nInfo.imat[id] = ((int)floorf(x/wdspp[1]) + 
			      (int)floorf(y/wdspp[1])*(int)wdspp[5] + 
			      (int)floorf(z/wdspp[1])*(int)(wdspp[5]*wdspp[5]));
#endif  

  mu = curand_uniform(&localState)*2-1;
  phi= curand_uniform(&localState)*2*PI;

  v[0] = sqrt(1.0-mu*mu)*cos(phi);
  v[1] = sqrt(1.0-mu*mu)*sin(phi);
  v[2] = mu; 

  int live=1;
#if defined(__CTALLY)
#endif
  //printf("[%2d],x=%.5f,pf=%.5f\n",id,DeviceMem.nInfo.pos_x[nid],pf);
  //for(istep=0;istep<devstep;istep++){
  while(live){
    l = -log(curand_uniform_double(&localState))*wdspp[2];
    //if(100>id) printf("  id=%d,l=%.5f\n",id,l);
    t = intersectbox(x,y,z,wdspp[0],wdspp[0],wdspp[0],v[0],v[1],v[2]);
    s = ((l/t+TEPSILON)<1)*l+((l/t+TEPSILON)>=1)*t;
    //s = ((l)<t)*l+((l)>=t)*t;
    //if(t<0) {printf("warning:t<0\n");                }
    //if(t>1.0e6) {printf("warning:t --> infinity \n");}
    x=x+s*v[0]; y=y+s*v[1]; z=z+s*v[2];
    live = (t>0)&&(t<1.0e6);//&&(x>0)&&(x<a)&&(y>0)&&(y<a)&&(z>0)&&(z<a);
    DeviceMem.nInfo.live[id] = live;
    s = live*s+(0==live)*t;
    if(  (t==s)||(x<=0)||(y<=0)||(z<=0)||(x>=wdspp[0])||(y>=wdspp[0])||(z>=wdspp[0])   ){//reflect
      //if slow, i can use the specific form for the box, which is changing sign of reflected component
      if((x-0)<TEPSILON){n[0]=-1.0;x = 0+TEPSILON;}
      if((y-0)<TEPSILON){n[1]=-1.0;y = 0+TEPSILON;}
      if((z-0)<TEPSILON){n[2]=-1.0;z = 0+TEPSILON;}
      if((wdspp[0]-x)<TEPSILON){n[0]= 1.0;x = wdspp[0]-TEPSILON;}
      if((wdspp[0]-y)<TEPSILON){n[1]= 1.0;y = wdspp[0]-TEPSILON;}
      if((wdspp[0]-z)<TEPSILON){n[2]= 1.0;z = wdspp[0]-TEPSILON;}
      add(v,n,-2*product(v,n));
      //printf("\n id=%2d, to (%.10e,%.10e,%.10e) along (%.10e,%.10e,%.10e)\n", blockDim.x * blockIdx.x + threadIdx.x,x,y,z,v[0],v[1],v[2]);
      //printf("id=%d, reflecting, time=%.3e\n",id,time);
      //fresh n
      n[0]=0; n[1]=0; n[2]=0;
    }  
    else{
#if defined(__CTALLY)
      nid = ((int)floorf(x/wdspp[1]) + (int)floorf(y/wdspp[1])*(int)wdspp[5] + (int)floorf(z/wdspp[1])*(int)(wdspp[5]*wdspp[5]));
      //DeviceMem.tally.cnt[ ((int)((int)(x/wdspp[1]) + (int)(y/wdspp[1])*wdspp[5] + (int) (z/wdspp[1])*wdspp[5]*wdspp[5]) )*gridDim.x*blockDim.x+id -shift ]+=1;
      DeviceMem.tally.cnt[ nid*gridDim.x*blockDim.x+id -shift ]+=1;
      //if(2==nid)
      //  printf("[id=%d]hitting 2 for %g times\n",id, DeviceMem.tally.cnt[ nid*gridDim.x*blockDim.x+id -shift ]);

#if defined(__CTALLY2)
      DeviceMem.cnt2_t[ nid*gridDim.x*blockDim.x+id -shift ]+=1;
#endif
      //if(34217==id) {printf("id=%d,ix=%.2f,iy=%.2f,iz=%.2f\n",id,floorf(x/wdspp[1]),floorf(y/wdspp[1]),floorf(z/wdspp[1]));}
      //if((nid<0)||(nid>=512)) printf("id=%d,index=%d,x=%.4f,y=%.4f,z=%.8e->%.8e,l=%.8e,t=%.8e,s=%.6e,%d,%d\n",id,nid,x,y,z-v[2]*s,z,l,t,s,(l<t),z<400.0);
#endif
      rnd = curand_uniform_double(&localState);
      if(rnd<(1-(wdspp[3]+wdspp[4]))){
	mu = curand_uniform(&localState)*2-1;
	phi= curand_uniform(&localState)*2*PI;
	v[0] = sqrt(1.0-mu*mu)*cos(phi);
	v[1] = sqrt(1.0-mu*mu)*sin(phi);
	v[2] = mu; 
	//if(100>id) printf("  id=%d, scatter to %.5f\n", id, mu);
      }
      else{
	live = 0;
	if(rnd>(1-wdspp[3])){ //fission
	  rnd = curand_uniform_double(&localState);
	  DeviceMem.nInfo.live[id] = 2*(rnd<=NU2)+3*(rnd>NU2);
	  //if(34217==id) printf("  id=%d, live[%d]= %d\n", id, id,DeviceMem.nInfo.live[id]);
	  //if(3<DeviceMem.nInfo.live[id]) printf("  id=%d, live[%d]= %d\n", id, id,DeviceMem.nInfo.live[id]);
          #if defined(__FTALLY2)
          DeviceMem.nInfo.imat[id-shift] = ((int)floorf(x/wdspp[1]) + 
                                      (int)floorf(y/wdspp[1])*(int)wdspp[5] + 
                                      (int)floorf(z/wdspp[1])*(int)(wdspp[5]*wdspp[5]));
          #endif  

	}
	else{  //rnd<Pc, capture, nothing to do
	  DeviceMem.nInfo.live[id] = 0;
	  //if(34217==id) printf("  id=%d, live[%d]= %d\n", id, id,DeviceMem.nInfo.live[id]);
	  //if(3<DeviceMem.nInfo.live[id]) printf("  id=%d, live[%d]= %d\n", id, id, DeviceMem.nInfo.live[id]);
	}
      }//end collision type

    }//end not leak
  }//end one history
  //}
#if defined(__CTALLY2)
  for(nid=0;nid<(wdspp[0]/wdspp[1]*wdspp[0]/wdspp[1]*wdspp[0]/wdspp[1]);nid++){
    DeviceMem.tally.cnt2[ nid*gridDim.x*blockDim.x+id -shift ] += (DeviceMem.cnt2_t[ nid*gridDim.x*blockDim.x+id -shift ]*DeviceMem.cnt2_t[ nid*gridDim.x*blockDim.x+id -shift ]);
  }
  //print contribution from each neutron at cell 2
  //if(DeviceMem.cnt2_t[ 2*gridDim.x*blockDim.x+id -shift ]){
  //  printf("%d\n", DeviceMem.cnt2_t[ 2*gridDim.x*blockDim.x+id -shift ]);
  //} 
#endif

  DeviceMem.nInfo.pos_x[id] = x;
  DeviceMem.nInfo.pos_y[id] = y;
  DeviceMem.nInfo.pos_z[id] = z;
#if defined(__MTALLY)
  DeviceMem.nInfo.imat[id] =  ((int)floorf(x/wdspp[1]) + 
			       (int)floorf(y/wdspp[1])*(int)wdspp[5] + 
			       (int)floorf(z/wdspp[1])*(int)(wdspp[5]*wdspp[5]) + 
			       DeviceMem.nInfo.imat[nid]*(int)(wdspp[5]*wdspp[5]*wdspp[5]) );
  //Note: wdspp[5] is num_bin in each dimension 
#endif

  DeviceMem.nInfo.rndState[id-shift] = localState; 
  //if(3<DeviceMem.nInfo.live[id]) printf("  id=%d, live[%d]= %d\n", id, id,DeviceMem);
  //if(34217==id) printf("  id=%d, live[%d]= %d, x=%.2f,y=%.2f,z=%.2f\n", id, id,DeviceMem.nInfo.live[id],x,y,z);
}
#endif //end if 3D

#if defined(__1D)
__device__ unsigned notleak(float x,float a){
  return (x>=0)&&(x<=a);
}
#if defined(__1D_VAC)
__global__ void history(MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize){
  float width = wdspp[0];
  float dx = wdspp[1];
  float mfp = wdspp[2];
  float Ps = 1-(wdspp[3]+wdspp[4]);
  float Pc = Ps+wdspp[4];
  float s;
  //try others when real simulation structure becomes clear
  //int idl = threadIdx.x;
  //id is the thread index
  //nid is the sampled index to get neutron position
  //in this scheme, normalization is realized by forcefully 
  //select gridsize neutrons from banksize neutrons
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  curandState localState = DeviceMem.nInfo.rndState[id];
  int nid = int(curand_uniform_double(&localState)*banksize);
  //extern __shared__ unsigned blockTerminated[];

  CMPTYPE rnd;
  float x = DeviceMem.nInfo.pos_x[nid];


  int dir = 1-2*int((curand_uniform_double(&localState))<=0.5);
  /* Copy state to local memory for efficiency */ 

  int newneu;
  unsigned live=1;
  //printf("[%2d],x=%.5f,pf=%.5f\n",id,DeviceMem.nInfo.pos_x[nid],pf);
  //for(istep=0;istep<devstep;istep++){
  while(live){
    s = -log(curand_uniform_double(&localState))*mfp;
    x = x+s*dir;

    live = notleak(x,width);
    //return true if not leak
    if(live){

      DeviceMem.tally.cnt[int(x/dx)*gridDim.x*blockDim.x+id]+=1;
    
      rnd = curand_uniform_double(&localState);
      if(rnd<Ps)
	dir = 1-2*int((curand_uniform_double(&localState))<=0.5);
      else{
	live = 0;
	if(rnd>Pc){ //fission
	  rnd = curand_uniform_double(&localState);
	  //newneu = 2*(rnd<=NU2)+3*(rand>NU2);
	  newneu = 1-2*(rnd<=NU2); //-1 --> 2 fission; +1 --> 3 fission
	  DeviceMem.nInfo.pos_y[id] = x*newneu;
	}
	else{  //rnd<Pc, capture, nothing to do
	  DeviceMem.nInfo.pos_y[id] = 0;
	}
      }//end collision type

    }//end not leak
  }//end one history
  //}
  //blockTerminated[idl] =1;// !live;
  
  /*Note: from now on, live does not indicate neutron but thread active */
  //blockActive[threadIdx.x] = (((terminated*2)*blockDim.x*gridDim.x + atomicAdd(Info.num_terminated_neutrons, terminated)) < num_src);
  //atomicAdd(Info.num_terminated_neutrons,!live);
  //Info.thread_active[id] =  blockDim.x*gridDim.x + *Info.num_terminated_neutrons < num_src;
  /* Copy state back to global memory */ 
  DeviceMem.nInfo.rndState[id] = localState; 

  /*
  else{
    blockTerminated[idl] = active;//0;
    //those old unlive neutrons must not be counted again
    //so, 0 instead of !live is used 
    //it was incorrect, above senario forgot to count leak neutron as terminated
  }

  //TODO: no need of such within block reduction for remaining()
  __syncthreads();
  live = blockDim.x>>1;
  while(live){
    if(idl<live)
      blockTerminated[idl] += blockTerminated[idl+live];
    __syncthreads();
    live>>=1;
  }
  if(0==idl){
    //reduction scheme depends on tally type
    //following is to count moderation times
    DeviceMem.block_terminated_neutrons[blockIdx.x] = blockTerminated[0];
  }
  */
}
#else
//_1D_Reflective
__global__ void history(MemStruct DeviceMem, unsigned num_src,int shift,unsigned banksize){
  float width = wdspp[0];
  float dx = wdspp[1];
  float mfp = wdspp[2];
  float Ps = 1-(wdspp[3]+wdspp[4]);
  float Pc = Ps+wdspp[4];
  float s, t, l;
  //try others when real simulation structure becomes clear
  //int idl = threadIdx.x;
  //id is the thread index
  //nid is the sampled index to get neutron position
  //in this scheme, normalization is realized by forcefully 
  //select gridsize neutrons from banksize neutrons
  int id = blockDim.x * blockIdx.x + threadIdx.x + shift;
  curandState localState = DeviceMem.nInfo.rndState[id-shift];
  int nid = int(curand_uniform_double(&localState)*banksize)+num_src;
  //extern __shared__ unsigned blockTerminated[];

  CMPTYPE rnd;
  float x = DeviceMem.nInfo.pos_x[nid];
#if defined(__FTALLY)
  DeviceMem.nInfo.imat[id] = int(floorf(x/dx));  
#endif  

  int dir = 1-2*int((curand_uniform_double(&localState))<=0.5);

  int newneu;
  unsigned live=1;
  //printf("[%2d],x=%.5f,pf=%.5f\n",id,DeviceMem.nInfo.pos_x[nid],pf);
  //for(istep=0;istep<devstep;istep++){
  while(live){
    l = -log(curand_uniform_double(&localState))*mfp;
    t = (width-x)*(1==dir)+(x)*((-1)==dir);
    s = ((l/t+TEPSILON)<1)*l+((l/t+TEPSILON)>=1)*t;
    x = x+s*dir;
    live = (t>0)&&(t<1.0e6);
    s = live*s+(0==live)*t;
    //while(!notleak(x,width)){
    //  x=((1==dir)*2*width+(-1==dir)*0-x);
    //  dir = 0-dir;
    //}
    if(  (t==s)||(x<=0)||(x>=wdspp[0]) ){//reflect
      if((x-0)<TEPSILON)       {dir= 1;x = 0       +TEPSILON;}
      if((wdspp[0]-x)<TEPSILON){dir=-1;x = wdspp[0]-TEPSILON;}
    }
    else{
#if defined(__CTALLY)
    DeviceMem.tally.cnt[ (int)(floorf(x/dx))*gridDim.x*blockDim.x+id -shift ]+=1;
#endif
    rnd = curand_uniform_double(&localState);
    if(rnd<Ps)
      dir = 1-2*int((curand_uniform_double(&localState))<=0.5);
    else{
      live = 0;
      if(rnd>Pc){ //fission
	rnd = curand_uniform_double(&localState);
	//newneu = 2*(rnd<=NU2)+3*(rand>NU2);
	newneu = 1-2*(rnd<=wdspp[6]); //-1 --> 2 fission; +1 --> 3 fission
	DeviceMem.nInfo.pos_y[id] = x*newneu;
      }
      else{  //rnd<Pc, capture, nothing to do
	DeviceMem.nInfo.pos_y[id] = 0;
      }
    }//end collision type
    }//end else reflect
  }//end one history
#if defined(__MTALLY)
  DeviceMem.nInfo.imat[id]= ( int(floorf(x/dx))  +DeviceMem.nInfo.imat[nid]*(int)(wdspp[5]));
  //Note: wdspp[5] is num_bin in each dimension 
#endif
  //}
  //blockTerminated[idl] =1;// !live;
  DeviceMem.nInfo.rndState[id-shift] = localState; 
}
#endif //end vac or reflective 1D
#endif //end if 1D

__global__ void reduce_sum_plus(int *threadcnt, int* cnt){
// reduce threadcnt[] to cnt[], cnt is updated by self increase
// this is used to count terminated neurtons
  /*reduce tally*/
  /*TODO:
    alternatives:
    1. only count for a block, saving global memory (acceess)
    2. count for each thread, saving time in thread wait
  */
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned idl = threadIdx.x;
  extern __shared__ int shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  int i;
  shared[idl] = threadcnt[id];
  __syncthreads();
  i = blockDim.x>>1;
  while(i){
    if(idl<i)
      shared[idl] += shared[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    //reduction scheme depends on tally type
    //following is to count moderation times
    cnt[blockIdx.x] += shared[0];
  }
}

__global__ void reduce_sum_equal(int* thread_active, int* active){
// reduce thread_active to active, active is updated without history
// this is used to count number of "live" threads
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned idl = threadIdx.x;
  extern __shared__ int shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  int i;
  shared[idl] = thread_active[id]; 
  __syncthreads();
  i = blockDim.x>>1;
  while(i){
    if(idl<i)
      shared[idl] += shared[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    active[blockIdx.x] = shared[0];
  }
}

__global__ void reduce_sum_equal(CMPTYPE* thread_active, CMPTYPE* active){
// reduce thread_active to active, active is updated without history
// this is used to count number of "live" threads
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned idl = threadIdx.x;
  extern __shared__ CMPTYPE sharedd[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  int i;
  sharedd[idl] = thread_active[id]; 
  __syncthreads();
  i = blockDim.x>>1;
  while(i){
    if(idl<i)
      sharedd[idl] += sharedd[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    active[blockIdx.x] = sharedd[0];
  }
}
