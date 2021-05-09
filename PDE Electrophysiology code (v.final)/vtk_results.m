%% Getting geometry from abaqus to export it to ensight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%% LOADING FILES OF NODES AND ELEM %%%%%%%%%%%%%%%%%%%%%%%%%

	% Loading nodes
load('mesh9.mat')
nodos = X;
nnode=length(X(:,1));
	% Loading elements
conectividades = R;
nelem=length(R(:,1));


%%%%%%%%%%%%%%%%%%%%%% WRITING HEADING FOR VTK FILE
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
        % printing heading to file
    f1=fopen('./results/nodes.vtk','w');
    f2=fopen('./results/elem.vtk','w');

    fprintf(f1,'# vtk DataFile Version 1.0\n');
    fprintf(f1,'ECM-CELL DIFFUSION-MECHANICS\n');
    fprintf(f1,'ASCII\n');
    fprintf(f1,'\n');
    fprintf(f1,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(f1,'%s %8i %s\n','POINTS', nnode ,'float');

    %%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%
    id = zeros(nnode, 1);
        % printing coordinates
    for i1=1:nnode
        fprintf(f1,'%14.8E %14.8E %14.8E\n',nodos(i1,1:3));
        id(i1)=i1;
    end
    fprintf(f1,'\n');

    %%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

    % printing connectivity

    fprintf(f2,'%s %8i %8i\n','CELLS', nelem , (nelem)*9);
    for i1=1:nelem
        new_conects(1,1)=conectividades(i1,1);
        new_conects(1,2)=conectividades(i1,2);
        new_conects(1,3)=conectividades(i1,3);
        new_conects(1,4)=conectividades(i1,4);
        new_conects(1,5)=conectividades(i1,5); 
        new_conects(1,6)=conectividades(i1,6); 
        new_conects(1,7)=conectividades(i1,7); 
        new_conects(1,8)=conectividades(i1,8); 
        fprintf(f2,'%4i  %10i  %10i  %10i %10i %10i  %10i  %10i %10i \n',8,new_conects(1,1)-1,...
            new_conects(1,2)-1,new_conects(1,3)-1,new_conects(1,4)-1, new_conects(1,5)-1,...
            new_conects(1,6)-1,new_conects(1,7)-1,new_conects(1,8)-1);
    end

    fprintf(f2,'\n');
    fprintf(f2,'%s %8i\n','CELL_TYPES', nelem);
    for i1=1:nelem
        fprintf(f2,' %4i ', 12);
    end
%%        
load('pde_electrophysiology.mat') 
 for i1=1:length(t)
    f3=fopen('./results/results.vtk','w');
    fprintf(f3,'%s %8i\n','POINT_DATA', nnode);
    fprintf(f3,'SCALARS GROUPS float\n');
    fprintf(f3,'LOOKUP_TABLE default\n');
    for i2=1:nnode
        fprintf(f3,'%8i\n', V_vtk(i2,i1) );
    end
    
    
%     fprintf(f3,'VECTORS orientation float\n');
%     for i2=1:nelem
%         fprintf(f,'%8i\n', [n1,n2,n3] );
%     end
    file_name1 = sprintf('cat ./results/nodes.vtk ./results/elem.vtk ./results/results.vtk > ./results/electro_%i.vtk',i1);
    system(file_name1);
 end