#include "optixmain.h"
char path_to_ptx[512];

void initialize_context(RTcontext context, int width, int n, int m, float *data, NeutronInfoStruct nInfo, multipole mp_para, CComplex<double>* wtable)
{
    /* Primary RTAPI objects */
    RTmaterial          material;
  
    /* Setup state */
#if defined(__HEXPRISM__)
    unsigned num_geobj = (1+3*n*(n+1))*(1+3*m*(m+1))+1 ;
    createContext( width,sqrt(3.f*m*m+3.f*m+1.f)*(n+1)*data[3]/*p*/,data[2]/*hh*/,num_geobj, context, nInfo,mp_para,wtable);
#else
    unsigned num_geobj = m*m*n*n*2+1 ;
    createContext( width,sqrt(2.0)*m*0.5*(n+2)*data[3]/*p*/,data[2]/*hh*/,num_geobj, context, nInfo,mp_para,wtable);
#endif

#if defined(__PRINTTRACK__)
    printf("%g,%g,%g,%g,%g\n",data[0],data[1],data[2],data[3],data[4]);
    printf("%d,%d,%d,%d,%d\n",width,n,m,0,0);
#endif
    createMaterial( context, &material);
    createInstances( context, material, data, n, m);

    /* Run */
    RT_CHECK_ERROR( rtContextValidate( context ) );
    RT_CHECK_ERROR( rtContextCompile( context ) );

    //disable optix time
    //printf("%g, %g, %g, %g, %g\n", time_elapsed, time_elapsed*1000.f/width/num_geobj, 0.f,0.f,0.f);
    // time cost (ms), average (us/geometry/ray)  
}

void createContext( int width, float R1, float Hh, unsigned num_geo, RTcontext context, NeutronInfoStruct nInfo, multipole mp_para, CComplex<double>* wtable)
{

    int id=0;
    RTprogram  ray_gen_program;
    RTprogram  miss_program;
    RTvariable only_one_ray_type;
    RTvariable epsilon;
    RTvariable max_depth;
    RTvariable var_R1, var_Hh, var_num; 

    RTvariable output_closest_buffer, output_current_buffer, output_live_buffer;
    RTvariable input_pos_x_buffer,input_pos_y_buffer,input_pos_z_buffer,
               input_dir_a_buffer,input_dir_p_buffer,
               input_id_buffer;
    RTbuffer   input_pos_x_buffer_obj, input_pos_y_buffer_obj, input_pos_z_buffer_obj,
               input_dir_a_buffer_obj, input_dir_p_buffer_obj,
               input_id_buffer_obj;
    RTbuffer            output_closest_buffer_obj;
    RTbuffer            output_current_buffer_obj;
    RTbuffer            output_live_buffer_obj;

    RTvariable dev_integers,
               mpdata, 
               wtable_buffer;
    RTbuffer   dev_integers_obj,
               mpdata_obj,
               wtable_buffer_obj;

    /* Setup context */
    RT_CHECK_ERROR( rtContextSetRayTypeCount( context, 2 ) );//TODO:type count /* shadow and radiance */
    RT_CHECK_ERROR( rtContextSetEntryPointCount( context, 1 ) );

    /*bind multipole parameters*/
    //It seems and makes sense that Optix does not care the size of the device memory, 
    //all multipole data are not available on host including array size, 
    //so, bufferSize for multipole parameters is meaningless
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "dev_integers", &dev_integers));
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &dev_integers_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( dev_integers_obj, RT_FORMAT_UNSIGNED_BYTE4)); 
    RT_CHECK_ERROR( rtBufferSetSize1D(dev_integers_obj, DEVINTS));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( dev_integers_obj, id, (CUdeviceptr)(mp_para.dev_integers)));
    RT_CHECK_ERROR( rtVariableSetObject( dev_integers, dev_integers_obj));

    RT_CHECK_ERROR( rtContextDeclareVariable( context, "mpdata", &mpdata));
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &mpdata_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( mpdata_obj,RT_FORMAT_USER )); 
    RT_CHECK_ERROR( rtBufferSetElementSize( mpdata_obj, sizeof(double)*2));
    RT_CHECK_ERROR( rtBufferSetSize1D(mpdata_obj, DEVINTS));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( mpdata_obj, id, (CUdeviceptr)(mp_para.mpdata)));
    RT_CHECK_ERROR( rtVariableSetObject( mpdata, mpdata_obj));

    RT_CHECK_ERROR( rtContextDeclareVariable( context, "wtable_buffer", &wtable_buffer));
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &wtable_buffer_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( wtable_buffer_obj,RT_FORMAT_USER )); 
    RT_CHECK_ERROR( rtBufferSetElementSize( wtable_buffer_obj, sizeof(double)*2));
    RT_CHECK_ERROR( rtBufferSetSize1D(wtable_buffer_obj, DEVINTS));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( wtable_buffer_obj, id, (CUdeviceptr)(wtable)));
    RT_CHECK_ERROR( rtVariableSetObject( wtable_buffer, wtable_buffer_obj));


    /*Declare variables*/
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "input_pos_x_buffer", &input_pos_x_buffer));
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "input_pos_y_buffer", &input_pos_y_buffer));
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "input_pos_z_buffer", &input_pos_z_buffer));
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "input_dir_p_buffer", &input_dir_p_buffer));
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "input_dir_a_buffer", &input_dir_a_buffer));

    RT_CHECK_ERROR( rtContextDeclareVariable( context, "input_id_buffer", &input_id_buffer));

    RT_CHECK_ERROR( rtContextDeclareVariable( context, "output_closest_buffer", &output_closest_buffer ) );
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "output_current_buffer", &output_current_buffer ) );
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "output_live_buffer", &output_live_buffer ) );

    RT_CHECK_ERROR( rtContextDeclareVariable( context, "max_depth", &max_depth ) );
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "only_one_ray_type", &only_one_ray_type ) );
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "scene_epsilon", &epsilon ) );
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "var_R1", &var_R1));
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "var_Hh", &var_Hh));
    RT_CHECK_ERROR( rtContextDeclareVariable( context, "var_num", &var_num));

    /*set variables: built-in types, rtVariableSetTYPE*/
    RT_CHECK_ERROR( rtVariableSet1i( max_depth, 10u ) );
    RT_CHECK_ERROR( rtVariableSet1ui( only_one_ray_type, 0u ) );
#if defined(__MANY__)
    RT_CHECK_ERROR( rtVariableSet1f( epsilon, 1.e-3f ) );
#else
    RT_CHECK_ERROR( rtVariableSet1f( epsilon, 1.e-3f ) );
#endif
    RT_CHECK_ERROR( rtVariableSet1f( var_R1, R1*0.9 ) );
    RT_CHECK_ERROR( rtVariableSet1f( var_Hh, Hh*0.9 ) );
    RT_CHECK_ERROR( rtVariableSet1ui( var_num, num_geo+1  ) );

    /* Render result buffer */
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &output_closest_buffer_obj) );
    RT_CHECK_ERROR( rtBufferSetFormat( output_closest_buffer_obj, RT_FORMAT_FLOAT ) );
    RT_CHECK_ERROR( rtBufferSetSize1D( output_closest_buffer_obj, width) );
    RT_CHECK_ERROR( rtBufferSetDevicePointer( output_closest_buffer_obj, id, (CUdeviceptr)(nInfo.d_closest)));
    RT_CHECK_ERROR( rtVariableSetObject( output_closest_buffer, output_closest_buffer_obj ) );

    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &output_current_buffer_obj) );
    RT_CHECK_ERROR( rtBufferSetFormat( output_current_buffer_obj, RT_FORMAT_UNSIGNED_BYTE4 ) );
    RT_CHECK_ERROR( rtBufferSetSize1D( output_current_buffer_obj, width) );
    RT_CHECK_ERROR( rtBufferSetDevicePointer( output_current_buffer_obj, id, (CUdeviceptr)(nInfo.imat)));
    //TODO: instances have not been assigned the materialID, instead cell ID is used
    RT_CHECK_ERROR( rtVariableSetObject( output_current_buffer, output_current_buffer_obj ) );

    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &output_live_buffer_obj) );
    RT_CHECK_ERROR( rtBufferSetFormat( output_live_buffer_obj, RT_FORMAT_UNSIGNED_BYTE4 ) );
    RT_CHECK_ERROR( rtBufferSetSize1D( output_live_buffer_obj, width) );
    RT_CHECK_ERROR( rtBufferSetDevicePointer( output_live_buffer_obj, id, (CUdeviceptr)(nInfo.live)));
    RT_CHECK_ERROR( rtVariableSetObject( output_live_buffer, output_live_buffer_obj ) );

    /* Input neutron id buffer*/
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &input_id_buffer_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( input_id_buffer_obj, RT_FORMAT_UNSIGNED_BYTE4)); 
    RT_CHECK_ERROR( rtBufferSetSize1D(input_id_buffer_obj, width));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( input_id_buffer_obj, id, (CUdeviceptr)(nInfo.id)));
    RT_CHECK_ERROR( rtVariableSetObject( input_id_buffer, input_id_buffer_obj));

    /* Input position buffer*/
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &input_pos_x_buffer_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( input_pos_x_buffer_obj, RT_FORMAT_FLOAT)); 
    RT_CHECK_ERROR( rtBufferSetSize1D(input_pos_x_buffer_obj, width));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( input_pos_x_buffer_obj, id, (CUdeviceptr)(nInfo.pos_x)));
    RT_CHECK_ERROR( rtVariableSetObject( input_pos_x_buffer, input_pos_x_buffer_obj));

    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &input_pos_y_buffer_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( input_pos_y_buffer_obj, RT_FORMAT_FLOAT)); 
    RT_CHECK_ERROR( rtBufferSetSize1D(input_pos_y_buffer_obj, width));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( input_pos_y_buffer_obj, id, (CUdeviceptr)(nInfo.pos_y)));
    RT_CHECK_ERROR( rtVariableSetObject( input_pos_y_buffer, input_pos_y_buffer_obj));

    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &input_pos_z_buffer_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( input_pos_z_buffer_obj, RT_FORMAT_FLOAT)); 
    RT_CHECK_ERROR( rtBufferSetSize1D(input_pos_z_buffer_obj, width));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( input_pos_z_buffer_obj, id, (CUdeviceptr)(nInfo.pos_z)));
    RT_CHECK_ERROR( rtVariableSetObject( input_pos_z_buffer, input_pos_z_buffer_obj));

    /* Input direction buffer*/
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &input_dir_a_buffer_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( input_dir_a_buffer_obj, RT_FORMAT_FLOAT)); 
    RT_CHECK_ERROR( rtBufferSetSize1D(input_dir_a_buffer_obj, width));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( input_dir_a_buffer_obj, id, (CUdeviceptr)(nInfo.dir_azimu)));
    RT_CHECK_ERROR( rtVariableSetObject( input_dir_a_buffer, input_dir_a_buffer_obj));
 
    RT_CHECK_ERROR( rtBufferCreateForCUDA( context, RT_BUFFER_INPUT, &input_dir_p_buffer_obj)); 
    RT_CHECK_ERROR( rtBufferSetFormat( input_dir_p_buffer_obj, RT_FORMAT_FLOAT)); 
    RT_CHECK_ERROR( rtBufferSetSize1D(input_dir_p_buffer_obj, width));
    RT_CHECK_ERROR( rtBufferSetDevicePointer( input_dir_p_buffer_obj, id, (CUdeviceptr)(nInfo.dir_polar)));
    RT_CHECK_ERROR( rtVariableSetObject( input_dir_p_buffer, input_dir_p_buffer_obj));

    /* Ray generation program */
#if defined(__MANY__)
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "pinhole_camera_many.ptx" );
#else
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "pinhole_camera_one.ptx" );
#endif
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "generate_ray", &ray_gen_program ) );
    //Declare and Set variables here if needed
    RT_CHECK_ERROR( rtContextSetRayGenerationProgram( context, 0, ray_gen_program ) );

    RTprogram  txs_program;
    RTvariable txs_program_var;
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "txs_one.ptx" );
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "xs_eval", &txs_program ) );
    RT_CHECK_ERROR( rtProgramDeclareVariable( ray_gen_program, "xs_eval", &txs_program_var));
    RT_CHECK_ERROR( rtVariableSetObject( txs_program_var, txs_program));


    /* Miss program */
#if defined(__MANY__)
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "constantbg_many.ptx" );
#else
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "constantbg_one.ptx" );
#endif
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "miss", &miss_program ) );
    RT_CHECK_ERROR( rtContextSetMissProgram( context, 0, miss_program ) );
}




void createGeometryBox( RTcontext context, RTgeometry* box, float* box_vertices )
{
	/*Geometry = boudingBox + intersection program*/
    RTprogram  box_intersection_program;
    RTprogram  box_bounding_box_program;
    RTvariable box_min_var;
    RTvariable box_max_var;

    RT_CHECK_ERROR( rtGeometryCreate( context, box ) );
    RT_CHECK_ERROR( rtGeometrySetPrimitiveCount( *box, 1u ) );

#if defined(__MANY__)
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "box_many.ptx" );
#else
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "box_one.ptx" );
#endif
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "box_bounds", &box_bounding_box_program ) );
    RT_CHECK_ERROR( rtGeometrySetBoundingBoxProgram( *box, box_bounding_box_program ) );
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "box_intersect", &box_intersection_program ) );
    RT_CHECK_ERROR( rtGeometrySetIntersectionProgram( *box, box_intersection_program ) );

	/*declare and set variable for "box_bounds" and "box_intersect" program*/
	/*boxmin and boxmax are used by both the programs*/

    RT_CHECK_ERROR( rtGeometryDeclareVariable( *box, "boxmin", &box_min_var ) );
    RT_CHECK_ERROR( rtGeometryDeclareVariable( *box, "boxmax", &box_max_var ) );
    RT_CHECK_ERROR( rtVariableSet3fv( box_min_var, box_vertices   ) );
    RT_CHECK_ERROR( rtVariableSet3fv( box_max_var, box_vertices+3 ) );
}

void createGeometryCylinder( RTcontext context, RTgeometry* cylinder,
                             float x1, float y1, float z1, float x2, float y2, float z2, float R){
	/*Geometry = boudingBox + intersection program*/
    RTprogram  cylinder_intersection_program;
    RTprogram  cylinder_bounding_box_program;
    RTvariable center1;
    RTvariable center2;
    RTvariable radius;

    RT_CHECK_ERROR( rtGeometryCreate( context, cylinder ) );
    RT_CHECK_ERROR( rtGeometrySetPrimitiveCount( *cylinder, 1u ) );

#if defined(__MANY__)
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "cyclinder_many.ptx" );
#else
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "cyclinder_one.ptx" );
#endif
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "cylinder_bounds", &cylinder_bounding_box_program ) );
    RT_CHECK_ERROR( rtGeometrySetBoundingBoxProgram( *cylinder, cylinder_bounding_box_program ) );
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "cylinder_intersect", &cylinder_intersection_program ) );
    RT_CHECK_ERROR( rtGeometrySetIntersectionProgram( *cylinder, cylinder_intersection_program ) );

    RT_CHECK_ERROR( rtGeometryDeclareVariable( *cylinder, "p1", &center1 ) );
    RT_CHECK_ERROR( rtGeometryDeclareVariable( *cylinder, "p2", &center2 ) );
    RT_CHECK_ERROR( rtGeometryDeclareVariable( *cylinder, "R", &radius ) );

    RT_CHECK_ERROR( rtVariableSet3f( center1, x1, y1, z1 ) );
    RT_CHECK_ERROR( rtVariableSet3f( center2, x2, y2, z2 ) );
    RT_CHECK_ERROR( rtVariableSet1f( radius,  R));
}

void createGeometrySphere( RTcontext context, RTgeometry* sphere, float* sphere_centeR )
{
	/*Geometry = boudingBox + intersection program*/
    RTprogram  sphere_intersection_program;
    RTprogram  sphere_bounding_box_program;
    RTvariable vsphere;

    RT_CHECK_ERROR( rtGeometryCreate( context, sphere ) );
    RT_CHECK_ERROR( rtGeometrySetPrimitiveCount( *sphere, 1u ) );

#if defined(__MANY__)
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "sphere_many.ptx" );
#else
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "sphere_one.ptx" );
#endif
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "sphere_bounds", &sphere_bounding_box_program ) );
    RT_CHECK_ERROR( rtGeometrySetBoundingBoxProgram( *sphere, sphere_bounding_box_program ) );
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "sphere_intersect", &sphere_intersection_program ) );
    RT_CHECK_ERROR( rtGeometrySetIntersectionProgram( *sphere, sphere_intersection_program ) );

    RT_CHECK_ERROR( rtGeometryDeclareVariable( *sphere, "sphere", &vsphere) );

    RT_CHECK_ERROR( rtVariableSet4fv( vsphere, sphere_centeR ) );
}



void createMaterial( RTcontext context, RTmaterial* material )
{
    /* Create our hit programs to be shared among all materials */

#if defined(__MANY__)
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "phong_many.ptx" );
    RT_CHECK_ERROR( rtMaterialCreate( context, material ) );
    RTprogram closest_hit_program;
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "closest_hit_radiance", &closest_hit_program ) );
    RT_CHECK_ERROR( rtMaterialSetClosestHitProgram( *material, 0, closest_hit_program ) );
#else
    sprintf( path_to_ptx, "%s/%s", "./obj/ptx", "phong_one.ptx" );
    RT_CHECK_ERROR( rtMaterialCreate( context, material ) );
    RTprogram any_hit_program;
    RT_CHECK_ERROR( rtProgramCreateFromPTXFile( context, path_to_ptx, "any_hit_shadow", &any_hit_program ) );
    RT_CHECK_ERROR( rtMaterialSetAnyHitProgram( *material, 0, any_hit_program ) );
#endif
}


void createInstances( RTcontext context, RTmaterial material, float *data, int n, int m)
{
    RTgeometry      cylinder;
    RTgroup         top_level_group;
    RTgeometrygroup geometrygroup;
    RTgeometryinstance instance;
    RTvariable      top_object;
    RTvariable      varID, matID;
    RTacceleration  acceleration,top_level_acceleration;
 
#if defined(__HEXPRISM__)
    int i = 0;
    int j = 0;
    int idir=0;
    int hdir=0;
    int index=0;
    float r1 = data[0];
    float r2 = data[1];
    float hh = data[2];
    float p  = data[3];
    float t  = data[4];
    float Hh = data[5];
    float phi = 0.f;
    float x0 = 0.f;
    float y0 = 0.f;
    float x = 0.f;
    float y = 0.f;

    /* create group to hold instances  */
    RT_CHECK_ERROR( rtGeometryGroupCreate( context, &geometrygroup ) );
    RT_CHECK_ERROR( rtGeometryGroupSetChildCount( geometrygroup, (1+3*n*(n+1))*(1+3*m*(m+1))+1 ) );
#if defined(__PRINTTRACK__)
    printf("%d, %d, %d, %d,%d\n",1+3*m*(m+1), 1+3*n*(n+1), (1+3*n*(n+1))*(1+3*m*(m+1))+1,0,0);
    // "constructing %d assemblies, each has %d fuel pins, which plus boundary sums to %d.
#endif  
    createGeometryCylinder(context, &cylinder, x0,y0,-hh, x0,y0,hh,r1);
    RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
    RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
    RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
    RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
    rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
    rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

    RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
    rtVariableSet1i(varID, ++index); 
    rtVariableSet1i(matID, 0); 

    //printf("cylinder: %d, x=%g, y=%g.\n", index,x0,y0);//pppp

    RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, 0, instance ) );

    for(i=1;i<=n;i++){
      x=p*i;
      y=0;
      createGeometryCylinder(context, &cylinder, x,y,-hh, x,y,hh,r1);
      RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
      RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
      RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
      RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
      rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
      rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

      RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
      rtVariableSet1i(varID, ++index); 
      rtVariableSet1i(matID, 0); 

      for(j=0;j<i*6-1;j++){
        idir=(j+i)/i-1;
        phi = (120+60*idir)*PI/180.f;
        x = x + p*cos(phi);
        y = y + p*sin(phi); 
        createGeometryCylinder(context, &cylinder, x,y,-hh, x,y,hh,r1);
        RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
        rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
        rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

        RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
        rtVariableSet1i(varID, ++index); 
        rtVariableSet1i(matID, 0); 

        ////printf("cylinder: %d, x=%g, y=%g.\n", index,x,y);//pppp
      } 
    }
    
    int k = 0;
    int h = 0;
    float l = sqrt(3.f)*(n+1)*p; 
    float phi2 = 0.f;
    float xp = 0.f; 
    float yp = 0.f; 

    for(k=1;k<=m;k++){
      xp=x0+l*k*cos(PI/6);
      yp=y0+l*k*sin(PI/6);
      createGeometryCylinder(context, &cylinder, xp,yp,-hh, xp,yp,hh,r1);
      RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
      RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
      RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
      RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
      rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
      rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

      RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
      rtVariableSet1i(varID, ++index); 
      rtVariableSet1i(matID, 0); 

      for(i=1;i<=n;i++){
        x=xp+p*i;
        y=yp+0;
        createGeometryCylinder(context, &cylinder, x,y,-hh, x,y,hh,r1);
        RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
        rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
        rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

        RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
        rtVariableSet1i(varID, ++index); 
        rtVariableSet1i(matID, 0); 
 
        for(j=0;j<i*6-1;j++){
          idir=(j+i)/i-1;
          phi = (120+60*idir)*PI/180.f;
          x = x + p*cos(phi);
          y = y + p*sin(phi); 
          createGeometryCylinder(context, &cylinder, x,y,-hh, x,y,hh,r1);
          RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
          RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
          RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
          RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
          rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
          rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

          RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
          rtVariableSet1i(varID, ++index); 
          rtVariableSet1i(matID, 0); 

          ////printf("cylinder: %d, x=%g, y=%g.\n", index,x,y);//pppp
        } 
      }
 
      for(h=0;h<6*k-1;h++){
        hdir=(h+k)/k-1;
        phi = (150+hdir*60)*PI/180.f;
        xp = xp + l*cos(phi);
        yp = yp + l*sin(phi);
        createGeometryCylinder(context, &cylinder, xp,yp,-hh, xp,yp,hh,r1);
        RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
        RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
        rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
        rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

        RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
        rtVariableSet1i(varID, ++index); 
        rtVariableSet1i(matID, 0); 

        ////printf("cylinder: %d, x=%g, y=%g.\n", index,xp,yp);//pppp
        for(i=1;i<=n;i++){
          x=xp+p*i;
          y=yp+0;
          createGeometryCylinder(context, &cylinder, x,y,-hh, x,y,hh,r1);
          RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
          RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
          RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
          RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
          rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
          rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

          RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
          rtVariableSet1i(varID, ++index); 
          rtVariableSet1i(matID, 0); 
   
          for(j=0;j<i*6-1;j++){
            idir=(j+i)/i-1;
            phi = (120+60*idir)*PI/180.f;
            x = x + p*cos(phi);
            y = y + p*sin(phi); 
            createGeometryCylinder(context, &cylinder, x,y,-hh, x,y,hh,r1);
            RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
            rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
            rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

            RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
            rtVariableSet1i(varID, ++index); 
            rtVariableSet1i(matID, 0); 

            ////printf("cylinder: %d, x=%g, y=%g.\n", index,x,y);//pppp
          } 
        }
   

      }    //end for h
    }      //end for k
    float R1 = sqrt(m*m+m+1.f/3.f)*l;
    createGeometryCylinder(context, &cylinder, x0,y0,-Hh, x0,y0,Hh,R1);
#else
    int i = 0;
    int k = 0;
    int ix=0,iy=0,kx=0,ky=0;
    int index=0;
    float r1 = data[0];
    float r2 = data[1];
    float hh = data[2];
    float p  = data[3];
    float t  = data[4];
    float Hh = data[5];
    float x0 = 0.f;
    float y0 = 0.f;
    float x = 0.f;
    float y = 0.f;

    /* create group to hold instances  */
    RT_CHECK_ERROR( rtGeometryGroupCreate( context, &geometrygroup ) );
    RT_CHECK_ERROR( rtGeometryGroupSetChildCount( geometrygroup, m*m*n*n*2+1 ));

#if defined(__PRINTTRACK__)
    printf("%d, %d, %d, %d,%d\n", m*m, n*n, m*m*n*n*2+1,0,0);
    // "constructing %d assemblies, each has %d fuel pins, which plus boundary sums to %d.
#endif
    float l = (n+2)*p; 
    float xp = 0.f; 
    float yp = 0.f; 

    for(k=0;k<m*m;k++){
        kx = k%m;
        ky = k/m;
        xp = x0 + l*kx;
        yp = y0 + l*ky;
        for(i=0;i<n*n;i++){
            ix = i%n;
            iy = i/n;
            x = xp + (1.5+ix)*p;
            y = yp + (1.5+iy)*p; 
            createGeometryCylinder(context, &cylinder, x,y,-hh+r2, x,y,hh-r2,r1);
            RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
            rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
            rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

            RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
            rtVariableSet1i(varID, ++index); 
            rtVariableSet1i(matID, 0); 

            createGeometryCylinder(context, &cylinder, x,y,-hh, x,y,hh,r2);
            RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
            RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
            rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
            rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

            RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
            rtVariableSet1i(varID, ++index); 
            rtVariableSet1i(matID, 0); 

        }  //end for i 
    }      //end for k
    float R1 = sqrt(2.)*m*l*0.5;
    createGeometryCylinder(context, &cylinder, x0+0.5*m*l,y0+0.5*m*l,-Hh, x0+0.5*m*l,y0+0.5*m*l,Hh,R1);
#endif
    RT_CHECK_ERROR( rtGeometryInstanceCreate( context, &instance ) );
    RT_CHECK_ERROR( rtGeometryInstanceSetGeometry( instance, cylinder ) );
    RT_CHECK_ERROR( rtGeometryInstanceSetMaterialCount( instance, 1 ) );
    RT_CHECK_ERROR( rtGeometryInstanceSetMaterial( instance, 0, material ) );
    rtGeometryInstanceDeclareVariable(instance, "geometryInstanceID", &varID);
    rtGeometryInstanceDeclareVariable(instance, "geometryMaterialID", &matID);

    RT_CHECK_ERROR( rtGeometryGroupSetChild( geometrygroup, index, instance ) );
    rtVariableSet1i(varID, ++index); 
    rtVariableSet1i(matID, 1); 

    /* create acceleration object for group and specify some build hints*/
    RT_CHECK_ERROR( rtAccelerationCreate(context,&acceleration) );
    RT_CHECK_ERROR( rtAccelerationSetBuilder(acceleration,BUILDER) );
    RT_CHECK_ERROR( rtAccelerationSetTraverser(acceleration,TRAVERSER) );
    RT_CHECK_ERROR( rtGeometryGroupSetAcceleration( geometrygroup, acceleration) );

    /* mark acceleration as dirty */
    RT_CHECK_ERROR( rtAccelerationMarkDirty( acceleration ) );


    RT_CHECK_ERROR( rtGroupCreate( context, &top_level_group ) );
    RT_CHECK_ERROR( rtGroupSetChildCount( top_level_group, 1 ) );
    RT_CHECK_ERROR( rtGroupSetChild( top_level_group, 0, geometrygroup ) );

    RT_CHECK_ERROR( rtContextDeclareVariable( context, "top_object", &top_object ) );
    RT_CHECK_ERROR( rtVariableSetObject( top_object, top_level_group ) );

    RT_CHECK_ERROR( rtAccelerationCreate( context, &top_level_acceleration ) );
    RT_CHECK_ERROR( rtAccelerationSetBuilder(top_level_acceleration,BUILDER) );
    RT_CHECK_ERROR( rtAccelerationSetTraverser(top_level_acceleration,TRAVERSER) );
    RT_CHECK_ERROR( rtGroupSetAcceleration( top_level_group, top_level_acceleration) );

    /* mark acceleration as dirty */
    RT_CHECK_ERROR( rtAccelerationMarkDirty( top_level_acceleration ) );
}


