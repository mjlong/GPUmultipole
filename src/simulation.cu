#include "simulation.h" 

extern __constant__ float wdspp[];
__global__ void initialize(MemStruct pInfo,float width, int banksize){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(id*7546861334684321478, id, id+14412078966483154, &(pInfo.nInfo.rndState[id]));

  neutron_sample(pInfo.nInfo, id,width);
  pInfo.nInfo.id[id] = id;
#if defined(__TALLY)
  pInfo.tally.cnt[id] = 0;
#endif 
  pInfo.nInfo.live[id] = 1*(id<banksize);
}

__device__ void neutron_sample(NeutronInfoStruct nInfo, unsigned id,float width){
  nInfo.live[id] = 1u;
  curandState state = nInfo.rndState[id];
  //TODO: source sampling should take settings dependent on geometry
#if defined(__1D)
  nInfo.pos_x[id] = width/PI*acos(1-2*curand_uniform_double(&state));//width*curand_uniform_double(&state);
#endif
#if defined(__3D)
#if defined(__TRAN)
  nInfo.dir_polar[id] = curand_uniform(&state)*2-1;
  nInfo.dir_azimu[id] = curand_uniform(&state)*PI*2;
  nInfo.d_closest[id] = 0.0; //used as time
  nInfo.pos_x[id] =width*curand_uniform_double(&state);
  nInfo.pos_y[id] =width*curand_uniform_double(&state); 
  nInfo.pos_z[id] =width*curand_uniform_double(&state); 
#else
  nInfo.pos_x[id] =width/PI*acos(1-2*curand_uniform_double(&state));//width*curand_uniform_double(&state);// 
  nInfo.pos_y[id] =width/PI*acos(1-2*curand_uniform_double(&state));//width*curand_uniform_double(&state);// 
  nInfo.pos_z[id] =width/PI*acos(1-2*curand_uniform_double(&state));//width*curand_uniform_double(&state);// 
#endif
#endif
  nInfo.rndState[id] = state;
#if defined(__WASTE)
  nInfo.energy[id] = STARTENE;
#endif
#if defined(__SCATTERPLOT)
  nInfo.energy[id] = nInfo.pos_z[id];
#endif
}

__device__ unsigned notleak(float x,float a){
  return (x>=0)&&(x<=a);
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

#if defined(__TRAN)
__global__ void history(MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned shift,float p2){
  float a = wdspp[0];
  float b=a;
  float c=a;

  float dx = wdspp[1];
  float mfp = wdspp[2];
  float Ps = 1-(wdspp[3]+wdspp[4]);
  float Pc = Ps+wdspp[4];

  float deltat = wdspp[2]/(wdspp[3]*2.5*wdspp[6]);
  //float deltat = 1.0/(wdspp[3]*2.5*wdspp[6]);
  float n[3] = {0.0, 0.0,0.0};
  float v[3];

  float l;//sampled path
  float t;//length to boundary
  float s;//min s
  float time=0.0;//real time

  int id = blockDim.x * blockIdx.x + threadIdx.x + shift;
  curandState localState = DeviceMem.nInfo.rndState[id];

  //extern __shared__ unsigned blockTerminated[];

  CMPTYPE rnd;
  int live    = DeviceMem.nInfo.live[id];
  //live can be -4,-3, -2,-1,0,1
  //            -4: over time more than two batches ago     --->-4
  //            -3: over time last batch                    --->-4
  //            -2: new fissioned                           ---> 1
  //            -1: terminated last time, still unlive      ---> 0
  //             0: terminated last time, still unlive
  //             1: continue
  float mu,phi; mu=0;phi=0;
  mu = DeviceMem.nInfo.dir_polar[id]*(1==live)+(curand_uniform(&localState)*2-1 )*(-2==live);
  phi= DeviceMem.nInfo.dir_azimu[id]*(1==live)+(curand_uniform(&localState)*2*PI)*(-2==live);

  live = ((-1!=live)&&(0!=live)&&(-3!=live)&&(-4!=live))+(-4)*(-3>=live); //-1 ---> 0 
  //live = (-2==live)||(1==live);  // now each neutron is the same

  if(1==live){
  float x = DeviceMem.nInfo.pos_x[id];//curand_uniform(&localState)*a;
  float y = DeviceMem.nInfo.pos_y[id];//curand_uniform(&localState)*b;
  float z = DeviceMem.nInfo.pos_z[id];//curand_uniform(&localState)*c;

  time    = DeviceMem.nInfo.d_closest[id]; 
  
  
  v[0] = sqrt(1.0-mu*mu)*cos(phi);
  v[1] = sqrt(1.0-mu*mu)*sin(phi);
  v[2] = mu; 


  while(1==live){
    l = -log(curand_uniform_double(&localState))*mfp;
  
    t = intersectbox(x,y,z,a,b,c,v[0],v[1],v[2]);
    //printf("\n id=%2d, from (%.10e,%.10e,%.10e) along (%.10e,%.10e,%.10e)\n", blockDim.x * blockIdx.x + threadIdx.x,x,y,z,v[0],v[1],v[2]);
    s = min(l,min(t,(deltat-time)*wdspp[6]));
    x=x+s*v[0]; y=y+s*v[1]; z=z+s*v[2]; time = time+s/wdspp[6];

    if(t==s){//reflect
      //if slow, i can use the specific form for the box, which is changing sign of reflected component
      if((x-0)<TEPSILON)
	n[0]=-1.0;
      if((y-0)<TEPSILON)
	n[1]=-1.0;
      if((z-0)<TEPSILON)
	n[2]=-1.0;

      if((a-x)<TEPSILON)
	n[0]= 1.0;
      if((b-y)<TEPSILON)
	n[1]= 1.0;
      if((c-z)<TEPSILON)
	n[2]= 1.0;
      add(v,n,-2*product(v,n));
      //printf("\n id=%2d, to (%.10e,%.10e,%.10e) along (%.10e,%.10e,%.10e)\n", blockDim.x * blockIdx.x + threadIdx.x,x,y,z,v[0],v[1],v[2]);
      //printf("id=%d, reflecting, time=%.3e\n",id,time);
      //fresh n
      n[0]=0; n[1]=0; n[2]=0;
    }  
    if(s==(deltat-time)*wdspp[6]){//time boundary
      live=-3;
      time = 0;
      //printf("id=%d, hitting time boundary, live=%d,time=%.3e\n",id,live,time);
    }

    //printf("[%2d],x=%.5f,pf=%.5f\n",id,DeviceMem.nInfo.pos_x[nid],pf);
    //for(istep=0;istep<devstep;istep++){

    if(s==l){//collison
#if defined(__TALLY)
      DeviceMem.tally.cnt[ (int(int(x/dx) + int(y/dx)*wdspp[5] + int (z/dx)*wdspp[5]*wdspp[5]) )*gridDim.x*blockDim.x+id]+=1;
#endif
    
      rnd = curand_uniform_double(&localState);
      if(rnd<Ps){
	v[2] = curand_uniform_double(&localState)*2-1;
	phi= curand_uniform_double(&localState)*2*PI;
	v[0] = sqrt(1-v[2]*v[2])*cos(phi);
	v[1] = sqrt(1-v[2]*v[2])*sin(phi);
	//printf("id=%d, scattering, live=%d,time=%.3e\n",id,live,time);
      }
      else{
	if(rnd>Pc){ //fission
	  rnd = curand_uniform_double(&localState);
	  live = 2*(rnd<=p2)+3*(rnd>p2);
	  //printf("id=%d,fission to %d, time=%.3e\n",id,live,time);
	}
	else{  //rnd<Pc, capture, nothing to do
	  live = -1;
	  //printf("id=%d, absorbed, live=%d,time=%.3e\n",id,live,time);
	}
      }//end collision type
    }//end if collision

  }//end one history
  //}
  //blockTerminated[idl] =1;// !live;
  
  /*Note: from now on, live does not indicate neutron but thread active */
  //blockActive[threadIdx.x] = (((terminated*2)*blockDim.x*gridDim.x + atomicAdd(Info.num_terminated_neutrons, terminated)) < num_src);
  //atomicAdd(Info.num_terminated_neutrons,!live);
  //Info.thread_active[id] =  blockDim.x*gridDim.x + *Info.num_terminated_neutrons < num_src;
  /* Copy state back to global memory */ 
  DeviceMem.nInfo.d_closest[id]= time;
  DeviceMem.nInfo.rndState[id] = localState; 
  DeviceMem.nInfo.pos_x[id] = x;
  DeviceMem.nInfo.pos_y[id] = y;
  DeviceMem.nInfo.pos_z[id] = z;
  DeviceMem.nInfo.dir_polar[id] = v[2];
  phi = v[0]/sqrt(1-v[2]*v[2]);//actually this is cosphi
  DeviceMem.nInfo.dir_azimu[id] = (v[1]>=0)*acos(phi)-(v[1]<0)*acos(phi);
  }//end if live
  //printf("id=%d, copying %d\n",id,live);
  DeviceMem.nInfo.live[id] = live;
  //printf("id=%d, %d copied \n",id,DeviceMem.nInfo.live[id]);
}
#else //3D steady State
__global__ void history(MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize){
  float a = wdspp[0];
  float b=a;
  float c=a;

  float dx = wdspp[1];
  float mfp = wdspp[2];
  float Ps = 1-(wdspp[3]+wdspp[4]);
  float Pc = Ps+wdspp[4];

  float n[3] = {0.0, 0.0,0.0};
  float v[3];

  float l;//sampled path
  float t;//length to boundary
  float s;

  int id = blockDim.x * blockIdx.x + threadIdx.x;
  curandState localState = DeviceMem.nInfo.rndState[id];
  int nid = int(curand_uniform_double(&localState)*banksize);
  //extern __shared__ unsigned blockTerminated[];

  CMPTYPE rnd;
  float x = DeviceMem.nInfo.pos_x[nid];
  float y = DeviceMem.nInfo.pos_y[nid];
  float z = DeviceMem.nInfo.pos_z[nid];
  float mu,phi; 
  mu = curand_uniform(&localState)*2-1;
  phi= curand_uniform(&localState)*2*PI;

  v[0] = sqrt(1.0-mu*mu)*cos(phi);
  v[1] = sqrt(1.0-mu*mu)*sin(phi);
  v[2] = mu; 

  unsigned live=1;
  //printf("[%2d],x=%.5f,pf=%.5f\n",id,DeviceMem.nInfo.pos_x[nid],pf);
  //for(istep=0;istep<devstep;istep++){
  while(live){
    l = -log(curand_uniform_double(&localState))*mfp;
    t = intersectbox(x,y,z,a,b,c,v[0],v[1],v[2]);
    s = (l<t)*l+(l>=t)*t;
    x=x+s*v[0]; y=y+s*v[1]; z=z+s*v[2];
    if(t==s){//reflect
      //if slow, i can use the specific form for the box, which is changing sign of reflected component
      if((x-0)<TEPSILON)
	n[0]=-1.0;
      if((y-0)<TEPSILON)
	n[1]=-1.0;
      if((z-0)<TEPSILON)
	n[2]=-1.0;

      if((a-x)<TEPSILON)
	n[0]= 1.0;
      if((b-y)<TEPSILON)
	n[1]= 1.0;
      if((c-z)<TEPSILON)
	n[2]= 1.0;
      add(v,n,-2*product(v,n));
      //printf("\n id=%2d, to (%.10e,%.10e,%.10e) along (%.10e,%.10e,%.10e)\n", blockDim.x * blockIdx.x + threadIdx.x,x,y,z,v[0],v[1],v[2]);
      //printf("id=%d, reflecting, time=%.3e\n",id,time);
      //fresh n
      n[0]=0; n[1]=0; n[2]=0;
    }  
    else{
#if defined(__TALLY)
      DeviceMem.tally.cnt[ (int(int(x/dx) + int(y/dx)*wdspp[5] + int (z/dx)*wdspp[5]*wdspp[5]) )*gridDim.x*blockDim.x+id  ]+=1;
#endif
      rnd = curand_uniform_double(&localState);
      if(rnd<Ps){
	mu = curand_uniform(&localState)*2-1;
	phi= curand_uniform(&localState)*2*PI;
	v[0] = sqrt(1.0-mu*mu)*cos(phi);
	v[1] = sqrt(1.0-mu*mu)*sin(phi);
	v[2] = mu; 
      }
      else{
	live = 0;
	if(rnd>Pc){ //fission
	  rnd = curand_uniform_double(&localState);
	  DeviceMem.nInfo.live[id] = 2*(rnd<=0.55)+3*(rnd>0.55);
	}
	else{  //rnd<Pc, capture, nothing to do
	  DeviceMem.nInfo.live[id] = 0;
	}
      }//end collision type

    }//end not leak
  }//end one history
  //}
  DeviceMem.nInfo.pos_x[id] = x;
  DeviceMem.nInfo.pos_y[id] = y;
  DeviceMem.nInfo.pos_z[id] = z;
  DeviceMem.nInfo.rndState[id] = localState; 
}
#endif //end tran or steady
#endif //end if 3D

#if defined(__1D)
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
	  //newneu = 2*(rnd<=0.55)+3*(rand>0.55);
	  newneu = 1-2*(rnd<=0.55); //-1 --> 2 fission; +1 --> 3 fission
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

    while(!notleak(x,width)){
      x=((1==dir)*2*width+(-1==dir)*0-x);
      dir = 0-dir;
    }
    DeviceMem.tally.cnt[int(x/dx)*gridDim.x*blockDim.x+id]+=1;
    
    rnd = curand_uniform_double(&localState);
    if(rnd<Ps)
      dir = 1-2*int((curand_uniform_double(&localState))<=0.5);
    else{
      live = 0;
      if(rnd>Pc){ //fission
	rnd = curand_uniform_double(&localState);
	//newneu = 2*(rnd<=0.55)+3*(rand>0.55);
	newneu = 1-2*(rnd<=0.55); //-1 --> 2 fission; +1 --> 3 fission
	DeviceMem.nInfo.pos_y[id] = x*newneu;
      }
      else{  //rnd<Pc, capture, nothing to do
	DeviceMem.nInfo.pos_y[id] = 0;
      }
    }//end collision type
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
