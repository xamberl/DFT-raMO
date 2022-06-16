function DFTraMOmovie()
%DFTraMOmovie: Generates the DFTraMOs for a series of desired target state using VASP results.
%
%    The program merges the DFT-raMO code of Vincent J. Yannello, with
%    the interface and target-generation functions of Daniel C.
%    Fredrickson's raMOmovie and some additional analysis functions
%    and debugs by Erdong Lu.
%
%    Copyright (C) 2020 Vincent J. Yannello, Erdong Lu and Daniel C. Fredrickson
%
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%



close all
fprintf('Welcome to the Fredrickson Group DFT raMO Movie Maker!\n\n');
stop = 0;

filename = 'sH_sc';
templatefile = 'template.xyz';
videofile = 'raMO.avi';
nruns = 0;
ntypes = 1;
nions = 1;
orblist_bytype = 1;
Emin = -17;
Emax = -5; 
Estep = 0.01;
broadening = 3;
DOSmax = 60;
E_Fermi = -7;
scriptcodes = [0 0 0 0 0 0 0];
scriptfiles{1,1} = '_blank_';
electron_counts = 0;
sphere_names{1,1} = 'XX';
sphere_radii(1,1) = 0.1;
sphere_colors(1,1) = 0; 
sphere_colors(1,2) = 0;
sphere_colors(1,3) = 0;
bond_names{1,1} = 'XX';
bond_names{1,2} = 'XX';
bond_maxlength(1,1) = 2.5;
bond_minlength(1,1) = 0.1;
bond_colors1(1,1) = 0; 
bond_colors1(1,2) = 0;
bond_colors1(1,3) = 0;
bond_colors2(1,1) = 0; 
bond_colors2(1,2) = 0;
bond_colors2(1,3) = 0;
bond_width(1,1) = 0.08;
npoly = 1;
poly_names{1,1} = 'XX';
poly_maxlength(1,1) = 2.5;
poly_minlength(1,1) = 0.1;
poly_colors1(1,1) = 0; 
poly_colors1(1,2) = 0;
poly_colors1(1,3) = 0;
poly_alpha(1,1) = 1;
poly_template{1,1} = 'template.xyz';

prompt = 'what is the SF?';
scalefactor = input(prompt)

while (stop == 0)    
    nspheres = size(sphere_radii);
    nspheres = nspheres(1);
    nbonds = size(bond_maxlength);
    nbonds = nbonds(1);
    npoly = size(poly_maxlength);
    npoly = npoly(1);
    fprintf('Menu options:    \n');
    fprintf('    1.  Change base file          [ %s ]\n',filename);      %  DONE
    fprintf('    2.  Change template file      [ %s ]\n',templatefile);  %  DONE
    fprintf('    3.  Change video file         [ %s ]\n',videofile);     %  DONE 
    fprintf('    4.  Add/edit raMO run         [ %d runs ]\n',nruns);    %  DONE
    fprintf('    5.  Change Nions by type      [ ');  %  DONE
    for j=1:ntypes
        fprintf('%d ',nions(1,j));
    end
    fprintf(']\n');
    fprintf('    6.  Change Norbitals by type  [ ');  %  DONE
    for j=1:ntypes
        fprintf('%d ',orblist_bytype(1,j));
    end
    fprintf(']\n');
    fprintf('    7.  Change Nelectrons by type [ ');  %  DONE
    for j=1:ntypes
        fprintf('%d ',electron_counts(1,j));
    end
    fprintf(']\n');
    fprintf('    8.  Change Fermi Energy       [ %f ]\n',E_Fermi);
    fprintf('    9.  Change Energy Range       [ %f %f %f %f %f ]\n',Emin,Emax,Estep,broadening,DOSmax);
    fprintf('   10.  Add/edit bonds            [ %d bond types]\n', nbonds);
    fprintf('   11.  Add/edit polyhedra        [ %d polyhedron types]\n', npoly);
    fprintf('   12.  Add/edit atomic spheres   [ %d atoms ]\n',nspheres);
    fprintf('   13.  Preview framework and DOS plots (Coming soon; for now use MullikenOptions.m)\n');
    fprintf('   14.  Load options\n');
    fprintf('   15.  Save options\n');
    fprintf('   16.  Run raMOmovie\n');
    fprintf('   17.  Quit\n\n');
    choice = input('Enter choice: ');
    fprintf('\n');

    if(choice==1) 
        filename = input('Enter new base file: ','s');
        fprintf('\n');
    end
    if(choice==2) 
        templatefile = input('Enter new template file: ','s');
        fprintf('\n');
    end
    if(choice==3) 
        videofile = input('Enter new video file: ','s');
        fprintf('\n');
    end
    if(choice==4) 
        fprintf('   Current list of raMO runs\n'); 
        nruns = size(scriptcodes);
        nruns = nruns(1);
        for n = 1:nruns
            fprintf('Run %d. ',n);
            if(scriptcodes(n,1) == 1)
                fprintf('Type:  Atomic Orbitals.  First site = %d.    Last site = %d.   First AO/atom:    %d.  Last AO/site:  %d.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5));    
            end
            if(scriptcodes(n,1) == 100)
                fprintf('Type:  Atomic Orbital bands.  First site = %d.    Last site = %d.   First AO/atom:    %d.  Last AO/site:  %d.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5));    
            end
            if(scriptcodes(n,1) == 2)
                fprintf('Type:  s-like cage.      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 200)
                fprintf('Type:  s-like cage bands.   First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 21)
                fprintf('Type:  s-like cage (p).      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 22)
                fprintf('Type:  s-like cage (sp3).      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 23)
                fprintf('Type:  s-like cage (sp3d2).      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 24)
                fprintf('Type:  interacting s-like cages.      First void = %d.    Last void = %d.     Void radius:    %f.  First MO = %d.    Last MO = %d.  Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5),scriptcodes(n,6),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 203)
                fprintf('Type:  s-like cage (sp3d2) bands.      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end

            if(scriptcodes(n,1) == 3)
                fprintf('Type:  Cluster MOs.      First center = %d.  Last center = %d.  Cluster radius:    %f.  First MO:    %d.  Last MO:  %d.  Centers file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5),scriptcodes(n,6),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 31)
                fprintf('Type:  Cluster MOs, Localized MOs.  First center = %d.  Last center = %d.  Cluster radius:    %f.  First MO:    %d.  Last MO:  %d.  Centers file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5),scriptcodes(n,6),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 41)
                fprintf('Type:  sp2d bonds.       First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 5)
                fprintf('Type:  Show remainders.  ');
            end
            if(scriptcodes(n,1) == 42)
                fprintf('Type:  sp2d remainders.  First atom = %d.    Last atom = %d.    ',scriptcodes(n,2), scriptcodes(n,3));    
            end
            if(scriptcodes(n,1) == 43)
                fprintf('Type:  sp2d (bis) bonds. First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 44)
                fprintf('Type:  sp2d (bis) remainders.  First atom = %d.    Last atom = %d.    ',scriptcodes(n,2), scriptcodes(n,3));    
            end
            if(scriptcodes(n,1) == 6)
                fprintf('Type:  Custom hybrids.      First atom = %d.    Last atom = %d.     Hybrids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 700)
                fprintf('Type:  Massive atomic Orbitals.  First site = %d.    Last site = %d.   First AO/atom:    %d.  Last AO/site:  %d.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5));
            end         
            fprintf(' Make movie = %d.\n',  scriptcodes(n,7));
        end
        nrun = input('Enter run to add/edit: ');
        if(nrun > 0)
          fprintf('Run types:  1 = Atomic orbital reconstructions. \n');
          fprintf('          100 = Atomic orbital reconstructions, with band creation. \n');
          fprintf('            2 = Reconstruction of s-like cage orbitals built from sp-hybrids on vertices. \n');
          fprintf('          200 = Reconstruction of s-like cage orbitals built from sp-hybrids on vertices, with filter and band creation. \n');
          fprintf('           21 = Reconstruction of s-like cage orbitals built from p orbitals on vertices. \n');
          fprintf('           22 = Reconstruction of s-like cage orbitals built from sp3-hybrids on vertices. \n');
          fprintf('           23 = Reconstruction of s-like cage orbitals built from sp3d2-hybrids on vertices. \n');
          fprintf('           24 = Reconstruction of MOs formed from s-like cage orbitals built from sp-hybrids on vertices. \n');
          fprintf('          203 = Reconstruction of s-like cage orbitals built from sp3d2-hybrids on vertices, with band creation. \n');
          fprintf('            3 = Reconstruction of MOs for cluster units defined by central points. \n');
          fprintf('           31 = Reconstruction of LCs of MOs for cluster units defined by central points. \n');
          fprintf('           41 = Reconstructions of bonds made with sp2d hydrids (on axis). \n');
          fprintf('           42 = Reconstructions of atomic orbitals remaining after sp2d hydridization (on axis). \n');
          fprintf('           43 = Reconstructions of bonds made with sp2d hydrids (bissecting axes). \n');
          fprintf('           44 = Reconstructions of atomic orbitals remaining after sp2d hydridization (bissecting axes). \n');
          fprintf('            5 = Show remainders. \n');
          fprintf('            6 = Custom hybrids. \n');
          fprintf('          700 = Massive atomic orbitals. \n');
          fprintf('\n');
          runtype = input('Enter run type: ');
          scriptcodes(nrun,1) = runtype;
          if(runtype == 1)
            scriptcodes(nrun,2) = input('Enter first atomic site to use: ');
            scriptcodes(nrun,3) = input('Enter  last atomic site to use: ');
            scriptcodes(nrun,4) = input('Enter first AO to use on each atom (1-9): ');
            scriptcodes(nrun,5) = input('Enter  last AO to use on each atom (1-9): ');
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = '_blank_';
          end
          if(runtype== 100)
            scriptcodes(nrun,2) = input('Enter first atomic site to use: ');
            scriptcodes(nrun,3) = input('Enter  last atomic site to use: ');
            scriptcodes(nrun,4) = input('Enter first AO to use on each atom (1-9): ');
            scriptcodes(nrun,5) = input('Enter  last AO to use on each atom (1-9): ');
            scriptcodes(nrun,6) = input('Scale factor for DOS of previous steps: ');
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = '_blank_';
          end

          if(runtype == 2)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
         if(runtype == 200)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = input('Scale factor for DOS of previous steps: ');
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
         end
         if(runtype == 22)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
          if(runtype == 23)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
          if(runtype == 24)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = input('Enter first MO to use: ');
            scriptcodes(nrun,6) = input('Enter last MO to use: ');;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
          if(runtype==203)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = input('Scale factor for DOS of previous steps: ');
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
          if(runtype == 21)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
          if(runtype == 6)
            scriptcodes(nrun,2) = input('Enter first atomic site to use: ');
            scriptcodes(nrun,3) = input('Enter  last atomic site to use: ');
            scriptcodes(nrun,4) = 0;
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with hydrid orbitals: ','s');
          end          
          if(runtype == 41)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
           if(runtype == 42)
            scriptcodes(nrun,2) = input('Enter first atom to use: ');
            scriptcodes(nrun,3) = input('Enter  last atom to use: ');
            scriptcodes(nrun,4) = 0;
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = '_blank_';
           end
          if(runtype == 43)
            scriptcodes(nrun,2) = input('Enter first void to use: ');
            scriptcodes(nrun,3) = input('Enter  last void to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with void centers: ','s');
          end
           if(runtype == 44)
            scriptcodes(nrun,2) = input('Enter first atom to use: ');
            scriptcodes(nrun,3) = input('Enter  last atom to use: ');
            scriptcodes(nrun,4) = 0;
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = '_blank_';
           end        
           if(runtype == 5)
            scriptcodes(nrun,2) = 0;
            scriptcodes(nrun,3) = 0;
            scriptcodes(nrun,4) = 0;
            scriptcodes(nrun,5) = 0;
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = '_blank_';
          end    
          if(runtype == 3)
            scriptcodes(nrun,2) = input('Enter first center to use: ');
            scriptcodes(nrun,3) = input('Enter  last center to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = input('Enter first MO to use on each cluster: ');
            scriptcodes(nrun,6) = input('Enter  last MO to use on each cluster: ');
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with cluster centers: ','s');
          end
          if(runtype == 31)
            scriptcodes(nrun,2) = input('Enter first center to use: ');
            scriptcodes(nrun,3) = input('Enter  last center to use: ');
            scriptcodes(nrun,4) = input('Enter  void radius: ');
            scriptcodes(nrun,5) = input('Enter first MO to use on each cluster: ');
            scriptcodes(nrun,6) = input('Enter  last MO to use on each cluster: ');
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = input('Enter file with cluster centers: ','s');
          end
          if(runtype == 700)
            fprintf('Do not record movie for this type of run! \n');
            scriptcodes(nrun,2) = input('Enter first atomic site to use: ');
            scriptcodes(nrun,3) = input('Enter  last atomic site to use: ');
            scriptcodes(nrun,4) = input('Enter first AO to use on each atom (1-9): ');
            scriptcodes(nrun,5) = input('Enter  last AO to use on each atom (1-9): ');
            scriptcodes(nrun,6) = 0;
            scriptcodes(nrun,7) = input('Record movie for this run [1 = yes]?  ');
            scriptfiles{nrun,1} = '_blank_';
          end

          fprintf('\nUpdated list of runs:  \n');
          nruns = size(scriptcodes);
          nruns = nruns(1);
          for n = 1:nruns
            fprintf('Run %d. ',n);
            if(scriptcodes(n,1) == 1)
                fprintf('Type:  Atomic Orbitals.  First site = %d.    Last site = %d.   First AO/atom:    %d.  Last AO/site:  %d.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5));    
            end
            if(scriptcodes(n,1) == 100)
                fprintf('Type:  Atomic Orbital bands.  First site = %d.    Last site = %d.   First AO/atom:    %d.  Last AO/site:  %d.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5));    
            end
            if(scriptcodes(n,1) == 2)
                fprintf('Type:  s-like cage.      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 200)
                fprintf('Type:  s-like cage bands.   First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 21)
                fprintf('Type:  s-like cage (p).      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 22)
                fprintf('Type:  s-like cage (sp3).      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 23)
                fprintf('Type:  s-like cage (sp3d2).      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 203)
                fprintf('Type:  s-like cage (sp3d2) bands.      First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 24)
                fprintf('Type:  interacting s-like cages.      First void = %d.    Last void = %d.     Void radius:    %f.  First MO = %d.    Last MO = %d.  Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5),scriptcodes(n,6),char(scriptfiles(n,1)));    
            end

            if(scriptcodes(n,1) == 3)
                fprintf('Type:  Cluster MOs.      First center = %d.  Last center = %d.  Cluster radius:    %f.  First MO:    %d.  Last MO:  %d.  Centers file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5),scriptcodes(n,6),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 31)
                fprintf('Type:  Cluster MOs, Localized MOs.  First center = %d.  Last center = %d.  Cluster radius:    %f.  First MO:    %d.  Last MO:  %d.  Centers file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5),scriptcodes(n,6),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 41)
                fprintf('Type:  sp2d bonds.       First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 42)
                fprintf('Type:  sp2d remainders.  First atom = %d.    Last atom = %d.   ',scriptcodes(n,2), scriptcodes(n,3));    
            end
            if(scriptcodes(n,1) == 43)
                fprintf('Type:  sp2d (bis) bonds. First void = %d.    Last void = %d.     Void radius:    %f.    Voids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),char(scriptfiles(n,1)));    
            end
            if(scriptcodes(n,1) == 44)
                fprintf('Type:  sp2d (bis) remainders.  First atom = %d.    Last atom = %d.    ',scriptcodes(n,2), scriptcodes(n,3));    
            end
            if(scriptcodes(n,1) == 5)
                fprintf('Type:  Show remainders.  ');
            end
            if(scriptcodes(n,1) == 6)
                fprintf('Type:  Custom hybrids.      First atom = %d.    Last atom = %d.     Hybrids file:  %s.  ',scriptcodes(n,2), scriptcodes(n,3),char(scriptfiles(n,1)));    
            end          
            if(scriptcodes(n,1) == 700)
                fprintf('Type:  Massive atomic Orbitals.  First site = %d.    Last site = %d.   First AO/atom:    %d.  Last AO/site:  %d.  ',scriptcodes(n,2), scriptcodes(n,3),scriptcodes(n,4),scriptcodes(n,5));
            end
            fprintf(' Make movie = %d.\n',  scriptcodes(n,7));
          end
        end
    end
    if(choice==5) 
        %fprintf('    5.  Change Nions by type      [ ');
        stop2 = 0;
        counter=0;
        while (stop2 == 0)
             fprintf('Enter number of atoms of type %d (0 if none)',counter+1);
             nion_temp = input(': ');
             if(nion_temp > 0)
                 counter = counter+1;
                 nions(1,counter) = nion_temp;
             else
                stop2 = 1;
                ntypes = counter;
             end
        end
        for j = 1:ntypes
            %   fprintf('    6.  Change Norbitals by type  [ ');
            %    for j=1:ntypes
            %        fprintf('%d ',orblist_bytype(1,j));
            %    end
            fprintf('Enter number of atomic orbitals for atoms of type %d',j);
            orblist_bytype(1,j)=input(': ');
        end
        for j = 1:ntypes
            %   fprintf('    6.  Change Norbitals by type  [ ');
            %    for j=1:ntypes
            %        fprintf('%d ',orblist_bytype(1,j));
            %    end
            fprintf('Enter number of electrons for atoms of type %d',j);
            electron_counts(1,j)=input(': ');
        end
    end
    if(choice==6) 
        for j = 1:ntypes
            %   fprintf('    6.  Change Norbitals by type  [ ');
            %    for j=1:ntypes
            %        fprintf('%d ',orblist_bytype(1,j));
            %    end
            fprintf('Enter number of atomic orbitals for atoms of type %d',j);
            orblist_bytype(1,j)=input(': ');
        end
    end
    if(choice==7) 
        for j = 1:ntypes
            %   fprintf('    6.  Change Norbitals by type  [ ');
            %    for j=1:ntypes
            %        fprintf('%d ',orblist_bytype(1,j));
            %    end
            fprintf('Enter number of electrons for atoms of type %d',j);
            electron_counts(1,j)=input(': ');
        end
    end
%    fprintf('    8.  Change Fermi Energy       [ %f ]\n',E_Fermi);
%    fprintf('    9.  Change DOS Parameters     [ %f %f %f %f %f ]\n',Emin,Emax,Estep,broadening,DOSmax);
    if(choice==8) 
           E_Fermi = input('Enter new Fermi Energy: ');
    end
    if(choice==9) 
           Emin = input('Enter new minimum E for DOS plot: ');
           Emax = input('Enter new maximum E for DOS plot: ');
           Estep = input('Enter new E step for DOS plot: ');
           broadening = input('Enter new Gaussian broadening parameter for DOS plot: ');
           DOSmax = input('Enter DOS maximum for plot: ');
    end
    if(choice==10)
        fprintf('\n');
        nbonds = size(bond_maxlength);
        nbonds = nbonds(1);
        for j = 1:nbonds
             fprintf('%d. %s [%f %f %f] to %s [%f %f %f].  Max length = %f. Min length = %f. Width = %f.\n',j,char(bond_names{j,1}),bond_colors1(j,1),bond_colors1(j,2),bond_colors1(j,3),char(bond_names{j,2}),bond_colors2(j,1),bond_colors2(j,2),bond_colors2(j,3),bond_maxlength(j,1),bond_minlength(j,1),bond_width(j,1)); 
        end
        fprintf('\n');
        nrun = input('Enter item to add/edit: ');
        if(nrun > 0)
            bond_names{nrun,1} = input('Enter atom name 1:  ','s');
            bond_colors1(nrun,1) = input('Enter color R code (0-255): ');
            bond_colors1(nrun,2) = input('Enter color G code (0-255): ');
            bond_colors1(nrun,3) = input('Enter color B code (0-255): ');
            bond_names{nrun,2} = input('Enter atom name 2:  ','s');
            if(strcmp(bond_names{nrun,1},bond_names{nrun,2})==1) 
                 bond_colors2(nrun,1)= bond_colors1(nrun,1);
                 bond_colors2(nrun,2)= bond_colors1(nrun,2);
                 bond_colors2(nrun,3)= bond_colors1(nrun,3);
            else     
                 bond_colors2(nrun,1) = input('Enter color R code (0-255): ');
                 bond_colors2(nrun,2) = input('Enter color G code (0-255): ');
                 bond_colors2(nrun,3) = input('Enter color B code (0-255): ');
            end
            bond_maxlength(nrun,1) = input('Enter max length: ');    
            bond_minlength(nrun,1) = input('Enter min length: ');    
            bond_width(nrun,1) = input('Enter bond width: ');
        end

        fprintf('\n');
        nbonds = size(bond_maxlength);
        nbonds = nbonds(1);
        for j = 1:nbonds
            fprintf('%d. %s [%f %f %f] to %s [%f %f %f].  Max length = %f. Min length = %f. Width = %f.\n',j,char(bond_names{j,1}),bond_colors1(j,1),bond_colors1(j,2),bond_colors1(j,3),char(bond_names{j,2}),bond_colors2(j,1),bond_colors2(j,2),bond_colors2(j,3),bond_maxlength(j,1),bond_minlength(j,1),bond_width(j,1)); 
        end
        fprintf('\n');
    end
   if(choice==11)
        fprintf('\n');
        npoly = size(poly_maxlength);
        npoly = npoly(1);
        for j = 1:npoly
             fprintf('%d. %s [%f %f %f].  Max length = %f. Min length = %f. Opacity = %f.  Template:  %s.\n',j,char(poly_names{j,1}),poly_colors1(j,1),poly_colors1(j,2),poly_colors1(j,3),poly_maxlength(j,1),poly_minlength(j,1),poly_alpha(j,1),char(poly_template{j,1})); 
        end
        fprintf('\n');
        nrun = input('Enter item to add/edit: ');
        if(nrun > 0)
            poly_names{nrun,1} = input('Enter central atom name:  ','s');
            poly_colors1(nrun,1) = input('Enter color R code (0-255): ');
            poly_colors1(nrun,2) = input('Enter color G code (0-255): ');
            poly_colors1(nrun,3) = input('Enter color B code (0-255): ');
            poly_maxlength(nrun,1) = input('Enter max length: ');    
            poly_minlength(nrun,1) = input('Enter min length: ');    
            poly_alpha(nrun,1) = input('Enter face opacity (0-1): ');
            poly_template{nrun,1} = input('Enter template file: ','s');
        end
        fprintf('\n');
        npoly = size(poly_maxlength);
        npoly = npoly(1);
        for j = 1:npoly
             fprintf('%d. %s [%f %f %f].  Max length = %f. Min length = %f. Opacity = %f.  Template:  %s.\n',j,char(poly_names{j,1}),poly_colors1(j,1),poly_colors1(j,2),poly_colors1(j,3),poly_maxlength(j,1),poly_minlength(j,1),poly_alpha(j,1),char(poly_template{j,1})); 
        end
        fprintf('\n');
    end
    if(choice==12)  % Change sphere radii
        fprintf('\n');
        nspheres = size(sphere_radii);
        nspheres = nspheres(1);
        for j = 1:nspheres
             fprintf('%d. Name = %s.  Radius = %f.  Color = [%f %f %f]\n',j,char(sphere_names{j,1}),sphere_radii(j,1),sphere_colors(j,1),sphere_colors(j,2),sphere_colors(j,3)); 
        end
        fprintf('\n');
        nrun = input('Enter item to add/edit: ');
        if(nrun > 0)
            sphere_names{nrun,1} = input('Enter atom name:  ','s');
            sphere_radii(nrun,1) = input('Enter sphere radius: ');
            sphere_colors(nrun,1) = input('Enter color R code (0-255): ');
            sphere_colors(nrun,2) = input('Enter color G code (0-255): ');
            sphere_colors(nrun,3) = input('Enter color B code (0-255): ');
        end
        fprintf('\n');
        nspheres = size(sphere_radii);
        nspheres = nspheres(1);
        for j = 1:nspheres
             fprintf('%d. Name = %s.  Radius = %f.  Color = [%f %f %f]\n',j,char(sphere_names{j,1}),sphere_radii(j,1),sphere_colors(j,1),sphere_colors(j,2),sphere_colors(j,3)); 
        end
        fprintf('\n');      
    end
    if(choice==14)  % Load options.
        fprintf('\n');
        ls
        fprintf('\n'); 
        options_file = input('Enter options file base name: ','s');
        num_options = strcat(options_file,'-num');
        string_options = strcat(options_file,'-str');
        sphere_options = strcat(options_file,'-sph');
        [sphere_names sphere_radii colorR colorG colorB] = textread(sphere_options,'%s %f %f %f %f');
        sphere_colors = [colorR colorG colorB];             
        bond_options = strcat(options_file,'-bnd');
        [bond_names1 colorR1 colorG1 colorB1 bond_names2 colorR2 colorG2 colorB2 bond_maxlength bond_minlength bond_width] = textread(bond_options,'%s %f %f %f %s %f %f %f %f %f %f');
        bond_colors1 = [colorR1 colorG1 colorB1];             
        bond_colors2 = [colorR2 colorG2 colorB2];
        bond_names = [bond_names1 bond_names2];

        poly_options = strcat(options_file,'-poly');
        [poly_names,poly_colors1R, poly_colors1G, poly_colors1B,poly_maxlength, poly_minlength, poly_alpha,poly_template] = textread(poly_options,'%s %f %f %f %f %f %f %s');
        poly_colors1 = [poly_colors1R, poly_colors1G, poly_colors1B];
        
        [option_name option_value] = textread(num_options,'%s %f');
        [str_option_name str_option_value] = textread(string_options,'%s %s');
        %  fprintf(f2,'basefile   %s\n',filename);
        %  fprintf(f2,'template   %s\n',templatefile);  
        %  fprintf(f2,'video      %s\n',videofile);    
        %  for n = 1:nruns
        %  fprintf(f2,'run_string %s\n',char(scriptfiles(n,1)));   
        %  end          
        nstr_options = size(str_option_value);
        nstr_options = nstr_options(1);
        filename = char(str_option_value(1,1));
        templatefile = char(str_option_value(2,1));
        videofile = char(str_option_value(3,1));
        scriptfiles = str_option_value(4:nstr_options,1);
        noptions = size(option_value);
        noptions = noptions(1);
        counter = 1;
        while counter < noptions+1
            if(strcmp(char(option_name(counter,1)),'Emin')==1)
                  Emin = option_value(counter,1);               
            end
            if(strcmp(char(option_name(counter,1)),'Emax')==1)
                  Emax = option_value(counter,1);               
            end
            if(strcmp(char(option_name(counter,1)),'Estep')==1)
                  Estep = option_value(counter,1);               
            end
            if(strcmp(char(option_name(counter,1)),'broadening')==1)
                  broadening = option_value(counter,1);               
            end
            if(strcmp(char(option_name(counter,1)),'DOSmax')==1)
                  DOSmax = option_value(counter,1);               
            end
            if(strcmp(char(option_name(counter,1)),'E_Fermi')==1)
                  E_Fermi = option_value(counter,1);               
            end
            if(strcmp(char(option_name(counter,1)),'ntypes')==1)
                  ntypes = option_value(counter,1);
                  for j2=1:ntypes
                        nions(1,j2) = option_value(counter+j2);
                  end
                  for j2=1:ntypes
                        orblist_bytype(1,j2) = option_value(counter+ntypes+j2);
                  end
                  for j2=1:ntypes
                        electron_counts(1,j2) = option_value(counter+2*ntypes+j2);
                  end
                  
                      
            end
            if(strcmp(char(option_name(counter,1)),'run')==1)
                   n = option_value(counter,1);
                   for j3 = 1:7  
                       scriptcodes(n,j3) = option_value(counter+j3,1);
                   end
                   nruns = n;
            end           
            counter = counter+1;
        end
    end
    
    if(choice==15)  % Save options.
        fprintf('\n');
        ls
        fprintf('\n'); 
        options_file = input('Enter options file base name: ','s');
        num_options = strcat(options_file,'-num');
        string_options = strcat(options_file,'-str');
        sphere_options = strcat(options_file,'-sph');
        f2 = fopen(sphere_options,'w');
        nspheres = size(sphere_radii);
        nspheres = nspheres(1);
        for j=1:nspheres
            fprintf(f2,'%s %f %f %f %f\n',char(sphere_names{j,1}),sphere_radii(j,1),sphere_colors(j,1), sphere_colors(j,2),sphere_colors(j,3));
        end
        fclose(f2);
        bond_options = strcat(options_file,'-bnd');
        f2 = fopen(bond_options,'w');
        nbonds = size(bond_maxlength);
        nbonds = nbonds(1);
        for j=1:nbonds
            fprintf(f2,'%s %f %f %f %s %f %f %f %f %f %f\n',char(bond_names{j,1}),bond_colors1(j,1), bond_colors1(j,2),bond_colors1(j,3),char(bond_names{j,2}),bond_colors2(j,1), bond_colors2(j,2),bond_colors2(j,3),bond_maxlength(j,1),bond_minlength(j,1), bond_width(j,1));
        end
        fclose(f2);
        poly_options = strcat(options_file,'-poly');
        f2 = fopen(poly_options,'w');
        npoly = size(poly_maxlength);
        npoly = npoly(1);
        for j=1:npoly
            fprintf(f2,'%s %f %f %f %f %f %f %s\n',char(poly_names{j,1}),poly_colors1(j,1), poly_colors1(j,2),poly_colors1(j,3),poly_maxlength(j,1),poly_minlength(j,1), poly_alpha(j,1),char(poly_template{j,1}));
        end
        fclose(f2);
        f2 = fopen(num_options,'w');
        fprintf(f2,'Emin        %f\n',Emin);
        fprintf(f2,'Emax        %f\n',Emax);
        fprintf(f2,'Estep       %f\n',Estep);
        fprintf(f2,'broadening  %f\n',broadening);
        fprintf(f2,'DOSmax      %f\n',DOSmax);
        fprintf(f2,'E_Fermi     %f\n',E_Fermi);
        fprintf(f2,'ntypes      %d\n',ntypes);
        for j=1:ntypes
          fprintf(f2,'nions       %d\n',nions(1,j)); 
        end
        for j=1:ntypes
          fprintf(f2,'orblist     %d\n',orblist_bytype(1,j)); 
        end
        for j=1:ntypes
          fprintf(f2,'nelectrons  %d\n',electron_counts(1,j)); 
        end
        nruns = size(scriptcodes);
        nruns = nruns(1);
        for n = 1:nruns
          fprintf(f2,'run         %d\n',n);
          fprintf(f2,'run_opt     %f\n',scriptcodes(n,1));
          fprintf(f2,'run_opt     %f\n',scriptcodes(n,2));
          fprintf(f2,'run_opt     %f\n',scriptcodes(n,3));
          fprintf(f2,'run_opt     %f\n',scriptcodes(n,4));
          fprintf(f2,'run_opt     %f\n',scriptcodes(n,5));
          fprintf(f2,'run_opt     %f\n',scriptcodes(n,6));
          fprintf(f2,'run_opt     %f\n',scriptcodes(n,7));
        end
        fclose(f2);
        f2 = fopen(string_options,'w');
        fprintf(f2,'basefile   %s\n',filename);
        fprintf(f2,'template   %s\n',templatefile);  
        fprintf(f2,'video      %s\n',videofile);    
        for n = 1:nruns
            fprintf(f2,'run_string %s\n',char(scriptfiles(n,1)));   
        end          
        fclose(f2);
    end
    
    if(choice==16)
        close all
        DOSparameters = [Emin Emax Estep broadening DOSmax];
        raMOmovie_run(filename, templatefile, videofile, nions, orblist_bytype, electron_counts, E_Fermi, DOSparameters,scriptcodes,scriptfiles,scalefactor); 
    end
    if(choice==17) 
        stop = 1; 
    end
end
end


function [new_psi_previous nelectrons_left]=reconstruct_targetsDFT(filename,templatefile,orblist,psi_target,movie_on, nelectrons_left, number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states,psi_previous,S_original,template_file,scalefactor)
      
%%%%%%%%    NOW CARRY OUT PROCESS ON DFT WAVEFUNCTIONS   %%%%%%%%%%%%%%%%
    %%% return variables: new_number_of_spin_up_states new_number_of_spin_down_states new_psi_previous
    n_targets= size(psi_target,2);
    nelectrons_left = nelectrons_left - 2*n_targets;
 output_name = strcat(filename,'-',int2str(nelectrons_left),'-psi_previous.mat');
% save(output_name,'new_psi_previous');   
if exist(output_name,'file')==2
    fprintf('   Found remainder file %s. Using these functions, and assuming that the functions of interest are coming later in the analysis.\n', output_name);
    load(output_name);
else
    atomsfile=strcat(filename,'-geo');
    [sH_name,sH_x,sH_y,sH_z]=textread(atomsfile,'%s %f %f %f');
    n_atoms_sc = size(orblist);
    n_atoms_sc = n_atoms_sc(2);
    cellfile=strcat(filename,'-cell');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    a_cellsc=a(1,:);
    b_cellsc=a(2,:);
    c_cellsc=a(3,:);
    
    use_fermi_energy=0;
    [atom_type]=outcar_reader; %Reads the OUTCAR file to determine identity of atoms of the cell
    number_of_atom_types=size(atom_type,2);
    atomic_numbers=atomlookup(atom_type); %Determines Atomic Numbers for each atom
    [acell,bcell,ccell,number_of_atoms_by_type,~,atom_position]=read_poscar(number_of_atom_types); %Reads in information from the POSCAR, atom locations, cell parameters.
    % acell, bcell, and ccell are column vectors
    astar=cross(bcell,ccell)/(dot(cross(bcell,ccell),acell));
    bstar=cross(ccell,acell)/(dot(cross(ccell,acell),bcell));
    cstar=cross(acell,bcell)/(dot(cross(acell,bcell),ccell));

    number_of_supercells_along_a = round(norm(a_cellsc)/norm(acell));
    number_of_supercells_along_b = round(norm(b_cellsc)/norm(bcell));
    number_of_supercells_along_c = round(norm(c_cellsc)/norm(ccell));


    %Calculating the overlap between each atomic orbital and the plane-waves
    currentorb=1;
    overlap_target_occupied=zeros(max(number_of_spin_up_states,number_of_spin_down_states),n_targets,number_of_spin_states);  
    [~,eht_symbol,~,eht_number_of_electrons,eht_number_of_orbitals,eht_N,eht_E,eht_Z1,eht_Z2,eht_n1,eht_n2]=readeht_parms; %Reads the Huckel parameters
    naos_sc = sum(orblist);
    H = zeros(naos_sc,naos_sc);
    for k2 = 1:n_atoms_sc
       frac_temp = linsolve([acell bcell ccell],[sH_x(k2,1) sH_y(k2,1) sH_z(k2,1)]');
       new_atom_positions(k2,1:3) = frac_temp';
       %  fprintf('Atom %d.  Cartesian:  %f %f %f.  Fractional:  %f %f %f.\n',k2,sH_x(k2,1), sH_y(k2,1), sH_z(k2,1), new_atom_positions(k2,1),new_atom_positions(k2,2),new_atom_positions(k2,3));
       [~,~,N,E,Z,n]=geteht_parms(sH_name(k2),eht_symbol, eht_number_of_orbitals, eht_number_of_electrons, eht_N, eht_E, eht_Z1, eht_Z2, eht_n1, eht_n2);
       prev_orb=sum(orblist(1:(k2-1)));
       if(orblist(k2) == 9) 
             H(prev_orb+1,prev_orb+1) = E(1,1);
             H(prev_orb+2,prev_orb+2) = E(2,1); 
             H(prev_orb+3,prev_orb+3) = E(2,1); 
             H(prev_orb+4,prev_orb+4) = E(2,1); 
             H(prev_orb+5,prev_orb+5) = E(3,1); 
             H(prev_orb+6,prev_orb+6) = E(3,1); 
             H(prev_orb+7,prev_orb+7) = E(3,1); 
             H(prev_orb+8,prev_orb+8) = E(3,1); 
             H(prev_orb+9,prev_orb+9) = E(3,1);
       end
       if(orblist(k2) == 4) 
             H(prev_orb+1,prev_orb+1) = E(1,1); 
             H(prev_orb+2,prev_orb+2) = E(2,1);
             H(prev_orb+3,prev_orb+3) = E(2,1);
             H(prev_orb+4,prev_orb+4) = E(2,1);
       end
       if(orblist(k2) == 1) 
             H(prev_orb+1,prev_orb+1) = E(1,1);
       end 
       % fprintf('Atom %d:  Hii_s = %f.  Hii_p = %f.  Hii_d = %f. \n',k2, E(1,1), E(2,1), E(3,1));
    end
    size(E)


    if number_of_spin_states==1
           spin_up_coefficients=1;
           spin_down_coefficients=1;
    else
        spin_up_coefficients=occupied_coefficients(:,1:number_of_spin_up_states);
        spin_down_coefficients=occupied_coefficients(:,number_of_spin_up_states+1:number_of_spin_up_states+number_of_spin_down_states);
    end
    for k=1:n_targets
      for k2=1:n_atoms_sc 
       prev_orb=sum(orblist(1:(k2-1)));
       if(norm(psi_target(prev_orb+1:prev_orb+orblist(k2),k))>0)
%         fprintf('Finding parameters for atom %d (%s).\n',char(sH_name(k2)));
%         sH_name
         [~,~,N,E,Z,n]=geteht_parms(sH_name(k2),eht_symbol, eht_number_of_orbitals, eht_number_of_electrons, eht_N, eht_E, eht_Z1, eht_Z2, eht_n1, eht_n2);
         Z=Z/0.52917721092; %Changing units from a.u. to angstroms.
         if number_of_spin_states==1
           spin_up_coefficients=1;
           spin_down_coefficients=1;
         end
         [overlap_target_temp,Emattemp]=calculate_target_occupied_overlap(kpoint_repeating,orblist(k2),number_of_planewaves,G,list_of_kpoints,new_atom_positions,astar,bstar,cstar,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states,number_of_occupied_states,occupied_coefficients,k2,N,E,Z,n,spin_up_coefficients,spin_down_coefficients);
         number_of_orbitals=size(overlap_target_temp,2);
         prev_orb=sum(orblist(1:(k2-1)));
         for ao_counter=1:orblist(k2)
             overlap_target_occupied(:,currentorb,:)=overlap_target_occupied(:,currentorb,:)+psi_target(prev_orb+ao_counter,k)*overlap_target_temp(:,ao_counter,:);
         end
       end
      end
      if k==1
           E_mat=psi_target(:,k)'*H*psi_target(:,k);
      else
           E_mat=[E_mat,psi_target(:,k)'*H*psi_target(:,k)];
      end
      currentorb=currentorb+1;     
    end
    %overlap_target_occupied=conj(overlap_target_occupied);
    %Generates the correct H and S matrices based on the nature of the analysis (number of spin states, number of remainder states).
    total_number_of_target_orbitals=currentorb-1;
    psi=zeros(size(overlap_target_occupied,1),size(overlap_target_occupied,1),number_of_spin_states);
    tr = size(psi_previous,1)-size(psi_previous,2); %tr = targets already constructed

    for j=1:number_of_spin_states
         if number_of_spin_states==2
            if j==1
                 H=psi_previous(1:number_of_spin_up_states,1:(number_of_spin_up_states-tr),j)'*overlap_target_occupied(1:number_of_spin_up_states,:,j)*diag(E_mat)*overlap_target_occupied(1:number_of_spin_up_states,:,j)'*psi_previous(1:number_of_spin_up_states,1:(number_of_spin_up_states-tr),j);
                 S=psi_previous(1:number_of_spin_up_states,1:(number_of_spin_up_states-tr),j)'*S_original(1:number_of_spin_up_states,1:number_of_spin_up_states)*psi_previous(1:number_of_spin_up_states,1:(number_of_spin_up_states-tr),j);
            else
            H=psi_previous(1:number_of_spin_down_states,1:(number_of_spin_down_states-tr),j)'*overlap_target_occupied(1:number_of_spin_down_states,:,j)*diag(E_mat)*overlap_target_occupied(1:number_of_spin_down_states,:,j)'*psi_previous(1:number_of_spin_down_states,1:(number_of_spin_down_states-tr),j);
            S=psi_previous(1:number_of_spin_down_states,1:(number_of_spin_down_states-tr),j)'*S_original(number_of_spin_up_states+1:end,number_of_spin_up_states+1:end)*psi_previous(1:number_of_spin_down_states,1:(number_of_spin_down_states-tr),j);
            end
            H=0.5*(H+H');
            S=0.5*(S+S');
            [psivect,entemp]=sort_eig(H,S+eps(eye(size(S,1))));
            %Sets up the new psi_previous or the final raMOs.
            if j==1
                tempup=psi_previous(1:number_of_spin_up_states,1:(number_of_spin_up_states-tr),j)*psivect;
                enup=entemp;
            else
                tempdown=psi_previous(1:number_of_spin_down_states,1:(number_of_spin_down_states-tr),j)*psivect;
                endown=entemp;
                if(number_of_spin_down_states-tr < 1)
                    tempdown = 0.0*tempup;
                    endown=0.0*enup;
                end
            end
         else
            H=psi_previous(:,:)'*overlap_target_occupied*diag(E_mat)*overlap_target_occupied'*psi_previous(:,:);
            S=psi_previous(:,:)'*S_original*psi_previous(:,:);
            H=0.5*(H+H');
            S=0.5*(S+S');
            [psivect,enup]=sort_eig(H,S+eps*eye(size(S,1)));
            %Sets up the new psi_previous or the final raMOs.
            tempup=psi_previous(:,:)*psivect;
         end
    end
    new_remainders = size(psi_previous,2)-n_targets;
    new_psi_previous=zeros(max(number_of_spin_up_states,number_of_spin_down_states),new_remainders,number_of_spin_states);
    if (number_of_spin_up_states-tr-n_targets) > 0
          new_psi_previous(1:number_of_spin_up_states,1:(number_of_spin_up_states-tr-n_targets),1)=tempup(:,n_targets+1:end);
    end
    fprintf('Number of spin down states remaining:  %d\n',number_of_spin_up_states-tr-n_targets);
    if number_of_spin_states==2
         if (number_of_spin_down_states-tr-n_targets) > 0
             new_psi_previous(1:number_of_spin_down_states,1:(number_of_spin_down_states-tr-n_targets),2)=tempdown(:,n_targets+1:end);
         end
         fprintf('Number of spin down states remaining:  %d\n',number_of_spin_down_states-tr-n_targets);
    end
    psiup=tempup(:,1:total_number_of_target_orbitals);
    if number_of_spin_states==2
        psidown=tempdown(:,1:total_number_of_target_orbitals);
    end
    %Plotting the raMO functions to .xsf files
    for j=1:number_of_spin_states
        if number_of_spin_states==2
            if j==1
                coeffprint=spin_up_coefficients;
                psiprint=psiup(1:number_of_spin_up_states,:);
                enprint=enup(1:total_number_of_target_orbitals);
                kpointsprint=list_of_kpoints(1:number_of_spin_up_states,:);
                kpointrepeatsprint=kpoint_repeating(1:number_of_spin_up_states);
            else
                coeffprint=spin_down_coefficients;
                psiprint=psidown(1:number_of_spin_down_states,:);
                enprint=endown(1:total_number_of_target_orbitals);
                kpointsprint=list_of_kpoints(number_of_spin_up_states+1:end,:);
                kpointrepeatsprint=kpoint_repeating(number_of_spin_up_states+1:end);
            end
        else
            coeffprint=occupied_coefficients;
            psiprint=psiup;
            enprint=enup(1:total_number_of_target_orbitals);
            kpointsprint=list_of_kpoints;
            kpointrepeatsprint=kpoint_repeating;
        end
        output_name = strcat(filename,'-',int2str(nelectrons_left),'-psi_previous.mat');
        save(output_name,'new_psi_previous');   
        if(movie_on > 0)   
            printing_size_limit=180000000;
            print_imaginary_results='off';
            degen_tol = 0.001;
            output_name = strcat(filename,'-',int2str(nelectrons_left),'-DFT');
            plot_wavefunction(G,psiprint,coeffprint,kpointsprint,atomic_numbers,enprint,kpointrepeatsprint,number_of_atom_types,number_of_supercells_along_a,number_of_supercells_along_b,number_of_supercells_along_c,printing_size_limit,output_name,print_imaginary_results,number_of_spin_states,j,degen_tol,template_file,scalefactor)
        end
    end
end
end



function plot_wavefunction(G,weights,occupied_coefficients,kpoint,atomic_numbers,en,kpoint_repeating,number_of_atom_types,number_of_supercells_along_a,number_of_supercells_along_b,number_of_supercells_along_c,printing_size_limit,output_name,print_imaginary_results,number_of_spin_states,spinstate,degen_tol,templatefile,scalefactor)
%Plots the wavefunctions to a .xsf file.
numprint=size(weights,2);
number_of_planewaves=size(G,1);
[nametemp,atomicno,xtemp,ytemp,ztemp]=textread(templatefile,'%s %d %f %f %f');


sizetemp=size(xtemp);
sizetemp=sizetemp(1);

[acell,bcell,ccell,atom_by_type,totnumatom,atom]=read_poscar(number_of_atom_types);
xtemp = xtemp; 
ytemp = ytemp;
ztemp = ztemp;

xyz_frac = zeros(sizetemp,3);
xyz_frac = [xtemp ytemp ztemp];
cellmatrix = [acell';bcell';ccell'];

for j = 1:sizetemp
   xyz_cart(j,1) =  xyz_frac(j,1)*cellmatrix(1,1)+xyz_frac(j,2)*cellmatrix(2,1)+xyz_frac(j,3)*cellmatrix(3,1);
   xyz_cart(j,2) =  xyz_frac(j,1)*cellmatrix(1,2)+xyz_frac(j,2)*cellmatrix(2,2)+xyz_frac(j,3)*cellmatrix(3,2);
   xyz_cart(j,3) =  xyz_frac(j,1)*cellmatrix(1,3)+xyz_frac(j,2)*cellmatrix(2,3)+xyz_frac(j,3)*cellmatrix(3,3);
end

x_min = min(xyz_frac(:,1));
x_max = max(xyz_frac(:,1));
y_min = min(xyz_frac(:,2));
y_max = max(xyz_frac(:,2));
z_min = min(xyz_frac(:,3));
z_max = max(xyz_frac(:,3));
xlength = (x_max-x_min)*norm(acell);
ylength = (y_max-y_min)*norm(bcell);
zlength = (z_max-z_min)*norm(ccell);

xrange = (x_max-x_min);
yrange = (y_max-y_min);
zrange = (z_max-z_min);
aratio=xlength/min([xlength ylength zlength]);
bratio=ylength/min([xlength ylength zlength]);
cratio=zlength/min([xlength ylength zlength]);

ratiofactor=(scalefactor*printing_size_limit/(aratio*bratio*cratio*number_of_planewaves))^(1/3);
ngfftsizea=floor(aratio*ratiofactor);
ngfftsizeb=floor(bratio*ratiofactor);
ngfftsizec=floor(cratio*ratiofactor);
orgshiftx = x_min*cellmatrix(1,1)+y_min*cellmatrix(2,1)+z_min*cellmatrix(3,1);
orgshifty = x_min*cellmatrix(1,2)+y_min*cellmatrix(2,2)+z_min*cellmatrix(3,2);
orgshiftz = x_min*cellmatrix(1,3)+y_min*cellmatrix(2,3)+z_min*cellmatrix(3,3);


fprintf('Voxel Grid size: %i %i %i.\n', ngfftsizea, ngfftsizeb, ngfftsizec);

%Sets up the otuput file names.
for i=1:numprint
  if number_of_spin_states==2
    if spinstate==1
      out{i}=strcat(output_name,'_',int2str(i),'up.xsf');
      f(i)=fopen(out{i}, 'w');
      if strcmp(print_imaginary_results,'on')==1
        out2{i}=strcat(output_name,'_',int2str(i),'upimag.xsf');
        f2(i)=fopen(out2{i}, 'w');
      end
    else
      out{i}=strcat(output_name,'_',int2str(i),'down.xsf');
      f(i)=fopen(out{i}, 'w');
      if strcmp(print_imaginary_results,'on')==1
        out2{i}=strcat(output_name,'_',int2str(i),'downimag.xsf');
        f2(i)=fopen(out2{i}, 'w');
      end
    end
  else
    out{i}=strcat(output_name,'_',int2str(i),'.xsf');
    f(i)=fopen(out{i}, 'w');
    if strcmp(print_imaginary_results,'on')==1
      out2{i}=strcat(output_name,'_',int2str(i),'imag.xsf');
      f2(i)=fopen(out2{i}, 'w');
    end
  end
  fprintf(f(i), ' DIM-GROUP\n 3 1\n PRIMVEC\n %f %f %f\n %f %f %f\n %f %f %f\n PRIMCOORD\n %i 1\n', xrange*acell(1), xrange*acell(2), xrange*acell(3), yrange*bcell(1), yrange*bcell(2), yrange*bcell(3), zrange*ccell(1), zrange*ccell(2), zrange*ccell(3), sizetemp);
  if strcmp(print_imaginary_results,'on')==1
    fprintf(f2(i), ' DIM-GROUP\n 3 1\n PRIMVEC\n %f %f %f\n %f %f %f\n %f %f %f\n PRIMCOORD\n %i 1\n', xrange*acell(1), xrange*acell(2), xrange*acell(3), yrange*bcell(1), yrange*bcell(2), yrange*bcell(3), zrange*ccell(1), zrange*ccell(2), zrange*ccell(3), sizetemp);
  end
end

%Generates the list of all atom positions, and the center for the output
count=1;
atomcenter=[0 0 0];

%Generates the grid of points
NX=zeros(ngfftsizea,ngfftsizeb,ngfftsizec);
NY=zeros(ngfftsizea,ngfftsizeb,ngfftsizec);
NZ=zeros(ngfftsizea,ngfftsizeb,ngfftsizec);
for i=1:ngfftsizea
  for j=1:ngfftsizeb
    for k=1:ngfftsizec
      NX(i,j,k)=xrange*(i-1)/(ngfftsizea-1)+x_min;
      NY(i,j,k)=yrange*(j-1)/(ngfftsizeb-1)+y_min;
      NZ(i,j,k)=zrange*(k-1)/(ngfftsizec-1)+z_min;
    end
  end
end

%Readjusts the cell parameters
currenttype=1;

%Prints the list of atoms
for i=1:sizetemp
  for l=1:numprint
    fprintf(f(l), '%i %f %f %f\n', atomicno(i), xyz_cart(i,1)-orgshiftx, xyz_cart(i,2)-orgshifty, xyz_cart(i,3)-orgshiftz);
    if strcmp(print_imaginary_results,'on')==1
      fprintf(f2(l), '%i %f %f %f\n', atomicno(i), xyz_cart(i,1)-orgshiftx, xyz_cart(i,2)-orgshifty, xyz_cart(i,3)-orgshiftz);
    end
  end
end
for l=1:numprint
  fprintf(f(l),'ATOMS\n');
  if strcmp(print_imaginary_results,'on')==1
    fprintf(f2(l),'ATOMS\n');
  end
end

numatomtemp=number_of_supercells_along_a*number_of_supercells_along_b*number_of_supercells_along_c*atom_by_type;
currenttype=1;

for i=1:sizetemp
  for l=1:numprint
    fprintf(f(l), '%i %f %f %f\n', atomicno(i), xyz_cart(i,1)-orgshiftx, xyz_cart(i,2)-orgshifty, xyz_cart(i,3)-orgshiftz);
    if strcmp(print_imaginary_results,'on')==1
      fprintf(f2(l), '%i %f %f %f\n', atomicno(i), xyz_cart(i,1)-orgshiftx, xyz_cart(i,2)-orgshifty, xyz_cart(i,3)-orgshiftz);
    end
  end
end
for l=1:numprint
  fprintf(f(l), 'BEGIN_BLOCK_DATAGRID3D\ndatagrids\nDATAGRID_3D_DENSITY\n %i %i %i\n 0.0 0.0 0.0\n %f %f %f\n %f %f %f\n %f %f %f\n', ngfftsizea, ngfftsizeb, ngfftsizec, xrange*acell(1), xrange*acell(2), xrange*acell(3), yrange*bcell(1), yrange*bcell(2), yrange*bcell(3), zrange*ccell(1), zrange*ccell(2), zrange*ccell(3));
  if strcmp(print_imaginary_results,'on')==1
    fprintf(f2(l), 'BEGIN_BLOCK_DATAGRID3D\ndatagrids\nDATAGRID_3D_DENSITY\n %i %i %i\n 0.0 0.0 0.0\n %f %f %f\n %f %f %f\n %f %f %f\n', ngfftsizea, ngfftsizeb, ngfftsizec, xrange*acell(1), xrange*acell(2), xrange*acell(3), yrange*bcell(1), yrange*bcell(2), yrange*bcell(3), zrange*ccell(1), zrange*ccell(2), zrange*ccell(3));
  end
end

%Generates the grid of wavefunctions
NX=reshape(NX,1,numel(NX));
NY=reshape(NY,1,numel(NY));
NZ=reshape(NZ,1,numel(NZ));
reducedkpoints=kpoint(kpoint_repeating==1,:);
kset=[];
count=1;
kpoint_repeating(size(kpoint_repeating,1)+1)=1;
for i=1:size(kpoint_repeating,1)-1
  kset=[kset,i];
  if kpoint_repeating(i+1)==1
    targetoverlap(:,:,count)=occupied_coefficients(:,kset)*weights(kset,:);
    kset=[];
    count=count+1;
  end
end
isosurf=zeros(ngfftsizea*ngfftsizeb*ngfftsizec,numprint);
size(targetoverlap)
for i=1:size(reducedkpoints,1)
% Attempting to save memory here. --DCF
  Gkx=bsxfun(@plus,G(:,1),reshape(reducedkpoints(i,1),1,size(reducedkpoints(i,1),1)));
  Gkx=bsxfun(@times,Gkx,NX);
  Gky=bsxfun(@plus,G(:,2),reshape(reducedkpoints(i,2),1,size(reducedkpoints(i,2),1)));
  Gkx=Gkx+bsxfun(@times,Gky,NY);
  clear Gky;
  Gkz=bsxfun(@plus,G(:,3),reshape(reducedkpoints(i,3),1,size(reducedkpoints(i,3),1)));
  Gkx=Gkx+bsxfun(@times,Gkz,NZ);
  clear Gkz;
%  Gkx=bsxfun(@plus,G(:,1),reshape(reducedkpoints(i,1),1,size(reducedkpoints(i,1),1)));
%  Gky=bsxfun(@plus,G(:,2),reshape(reducedkpoints(i,2),1,size(reducedkpoints(i,2),1)));
%  Gkz=bsxfun(@plus,G(:,3),reshape(reducedkpoints(i,3),1,size(reducedkpoints(i,3),1)));
%  Gkx=bsxfun(@times,Gkx,NX);
%  Gky=bsxfun(@times,Gky,NY);
%  Gkz=bsxfun(@times,Gkz,NZ);
%  planewave=exp(pi*2i*(Gkx+Gky+Gkz));
  planewave=exp(pi*2i*(Gkx));
  isosurf=isosurf+planewave.'*targetoverlap(:,:,i);
  fprintf('K-point Number %i Printed.\n', i);
  clear planewave;
end
clear Gkx;

isosurf=reshape(isosurf,ngfftsizea,ngfftsizeb,ngfftsizec,numprint);

%The Wavefunctions read in are in general complex. Now we need to rotate them back onto the real axis. This code does this, as well as deals with degenerate states 
E=en;
E(numprint+1)=10000;
numdegen=0;
for j=1:numprint
  if numdegen>1
    numdegen=numdegen-1;
    continue;
  end
  numdegen=0;
  temp=0;
  while numdegen==temp
    if abs(E(j+temp+1)-E(j))<degen_tol*abs(E(j))
      temp=temp+1;
      numdegen=numdegen+1;
    else
      numdegen=numdegen+1;
    end
  end
  tempmatrix=reshape(isosurf(:,:,:,j:j+numdegen-1),ngfftsizea*ngfftsizeb*ngfftsizec,numdegen);
  randos=randi(ngfftsizea*ngfftsizeb*ngfftsizec,1000,1);
  realimag=[real(tempmatrix(randos,:)),imag(tempmatrix(randos,:))];
  [combo,E_combo]=eig(realimag'*realimag);
  [~,v]=sort(diag(abs(real(E_combo))));
  combo=combo(:,v(1:numdegen));
  combo2=complex(combo(numdegen+1:2*numdegen,:),combo(1:numdegen,:));
  tempmatrix=tempmatrix*combo2;
  isosurf2(:,:,:,j:j+numdegen-1)=reshape(tempmatrix,ngfftsizea,ngfftsizeb,ngfftsizec,numdegen);
end
for j=1:numprint
  fprintf(f(j), '%f\n', isosurf2(:,:,:,j));
  if strcmp(print_imaginary_results,'on')==1
    fprintf(f2(j), '%f\n', imag(isosurf2(:,:,:,j)));
  end
end
for j=1:numprint
  fprintf(f(j), 'END_DATAGRID_3D\n END_BLOCK_DATAGRID3D');
  if strcmp(print_imaginary_results,'on')==1
    fprintf(f2(j), 'END_DATAGRID_3D\n END_BLOCK_DATAGRID3D');
  end
end
end


function raMOmovie_run(filename, templatefile, videofile, nions, orblist_bytype, electron_counts, E_Fermi, DOSparameters,scriptcodes,scriptfiles,scalefactor)
Eshift = E_Fermi;
% E_Fermi = -7.712518;
% DOSparameters = [Emin,Emax,Estep,broadening,DOSmax]
%
%
%
%
%
%  scriptcodes:  gives set of targets to analyze in order of type
%  scriptfiles:  filenames needed for each run
%
%
%
%
%

n_scriptruns = size(scriptcodes);
n_scriptruns = n_scriptruns(1);

  cellfile=strcat(filename,'-cell');
  a=zeros(3,3);
  [a1,a2,a3]=textread(cellfile,' %f %f %f');
  a=[a1,a2,a3];
  a_cellsc=a(1,:);
  b_cellsc=a(2,:);
  c_cellsc=a(3,:);
  
%%%%  INSERTING INITIALIZION OF DFT-raMO VARIABLES HERE  %%%%

% VARIABLES IN DFT-raMO CODE TO BE DEFINED.
%number_of_supercells_along_a
%number_of_supercells_along_b
%number_of_supercells_along_c
%occupied_filling_level
%use_fermi_energy
%fermi_energy
%printing_size_limit
%print_imaginary_results
%degen_tol
Emin=        DOSparameters(1,1);
Emax=        DOSparameters(1,2);
Estep=       DOSparameters(1,3);
broadening=  DOSparameters(1,4);
DOSmax =     DOSparameters(1,5);

use_fermi_energy=0;

[atom_type]=outcar_reader; %Reads the OUTCAR file to determine identity of atoms of the cell
number_of_atom_types=size(atom_type,2);
[acell,bcell,ccell,number_of_atoms_by_type,~,atom_position]=read_poscar(number_of_atom_types); %Reads in information from the POSCAR, atom locations, cell parameters.


%  DETERMINE SUPERCELL DIMENSIONS FROM COMPARISON OF sH SUPERCELL AND VASP
%  CELL
number_of_supercells_along_a = round(norm(a_cellsc)/norm(acell));
number_of_supercells_along_b = round(norm(b_cellsc)/norm(bcell));
number_of_supercells_along_c = round(norm(c_cellsc)/norm(ccell));

fprintf('Comparing cells from sH and VASP:  sH geometry represents %d x %d x %d supercell of VASP cell.\n',number_of_supercells_along_a,number_of_supercells_along_b,number_of_supercells_along_c) 

atomic_numbers=atomlookup(atom_type); %Determines Atomic Numbers for each atom
[~,eht_symbol,~,eht_number_of_electrons,eht_number_of_orbitals,eht_N,eht_E,eht_Z1,eht_Z2,eht_n1,eht_n2]=readeht_parms; %Reads the Huckel parameters

%Reads the readableparams.txt file, or makes it and then reads it
input_name = 'GCOEFF.txt';
occupied_filling_level = 0.5;
name_of_readableparams=strcat('readableparams_',num2str(Emin),'_',num2str(Emax),'.txt');
if exist(name_of_readableparams,'file')==2
  [number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states]=read_readableparams(name_of_readableparams);
else
  make_readableparams(input_name,name_of_readableparams,occupied_filling_level,use_fermi_energy,Emin,Emax);
  [number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states]=read_readableparams(name_of_readableparams);
end

%%%  NOW HAVE VASP WAVEFUNCTIONS IN MEMORY  %%%

%Separates out the two spin states or if there is only one spin state, sets a variable for convenience sake.
if number_of_spin_states==1;
  number_of_spin_up_states=number_of_occupied_states;
else
  spin_up_coefficients=occupied_coefficients(:,1:number_of_spin_up_states);
  spin_down_coefficients=occupied_coefficients(:,number_of_spin_up_states+1:number_of_spin_up_states+number_of_spin_down_states);
end

%Generates the Overlap Matrix for the occupiedstates
S=make_overlap_mat(occupied_coefficients,kpoint_repeating);
S=0.5*(S+S');
S_original=S;

%Calculating astar, bstar, and cstar for later.
astar=cross(bcell,ccell)/(dot(cross(bcell,ccell),acell));
bstar=cross(ccell,acell)/(dot(cross(ccell,acell),bcell));
cstar=cross(acell,bcell)/(dot(cross(acell,bcell),ccell));

%Reads in the kpoint information from the KPOINTS file.
max_number_of_kpoints=read_kpoints;
previous_target_number=0; %Sets the previous number of target states to 0.

%psi_previous = columns are current wavefunctions with the elements along
%down the column being the coefficients from the original occupied states.
%To start out, this is just the identity matrix.
%This what will be updated at the end of each reconstruct_target run.

psi_previous(:,:,1)=eye(max(number_of_spin_up_states,number_of_spin_down_states),max(number_of_spin_up_states,number_of_spin_down_states));
if number_of_spin_states==2
   psi_previous(:,:,2)=eye(max(number_of_spin_up_states,number_of_spin_down_states),max(number_of_spin_up_states,number_of_spin_down_states));
end

%%%%%%%%%%%%%%%%%%END OF CODE FROM DFT-raMO  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ntypes=size(nions);
ntypes=ntypes(2);

orblistW = [orblist_bytype(1,1)*ones(1,nions(1,1))];

for j = 2:ntypes
    orblistW = [orblistW orblist_bytype(1,j)*ones(1,nions(1,j))];
end


norbs = sum(orblistW);

nelectrons=0;
for j = 1:ntypes
    nelectrons = nelectrons+ nions(1,j)*electron_counts(1,j);
end


cellfile=strcat(filename,'-cell');
geofile=strcat(filename,'-geo');


orb_counter=0;

fprintf('Number of electrons from compound formula:     %d\n',nelectrons);
nelectrons_left = nelectrons;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for script_run = 1:n_scriptruns
    fprintf('STARTING RUN %d.\n',script_run);
%  FOCUS ON ATOMIC ORBITAL RECONSTRUCTIONS AT FIRST.
    if (scriptcodes(script_run,1) == 1) 
        first_site  = scriptcodes(script_run,2);
        second_site = scriptcodes(script_run,3);
        first_orbital = scriptcodes(script_run,4);
        last_orbital = scriptcodes(script_run,5);
        movie_on = scriptcodes(script_run,7);
        for siteno =   first_site:second_site
           psi_target = make_target_AO(norbs, siteno, orblistW, first_orbital:last_orbital);
           fprintf('   RUN %d:  Targets created for site %d orbitals %d:%d.\n',script_run, siteno,first_orbital,last_orbital);

           [psi_previous nelectrons_left]=reconstruct_targetsDFT(filename,templatefile,orblistW, psi_target,movie_on, nelectrons_left,number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states, psi_previous,S_original,templatefile,scalefactor);
        end 
    end

    if (scriptcodes(script_run,1) == 700)
        first_site  = scriptcodes(script_run,2);
        second_site = scriptcodes(script_run,3);
        first_orbital = scriptcodes(script_run,4);
        last_orbital = scriptcodes(script_run,5);
        movie_on = scriptcodes(script_run,7);
        psi_target = make_target_massAO(norbs, first_site, second_site, orblistW, first_orbital:last_orbital);
        fprintf('   RUN %d:  Targets created for atoms %d:%d orbitals %d:%d.\n',script_run, first_site, second_site, first_orbital, last_orbital);
        [psi_previous nelectrons_left]=reconstruct_targetsDFT(filename,templatefile,orblistW, psi_target,movie_on, nelectrons_left,number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states, psi_previous,S_original,templatefile,scalefactor);
    end

%  CODE FOR ANALYZING OTHER TYPES OF TARGETS TO BE MODIFIED LATER 
%     if (scriptcodes(script_run,1) == 100)
%         first_site  = scriptcodes(script_run,2);
%         second_site = scriptcodes(script_run,3);
%         first_orbital = scriptcodes(script_run,4);
%         last_orbital = scriptcodes(script_run,5);
%         DOS_scale = scriptcodes(script_run,6);
%         movie_on = 0;
%         nraMOs=0;
%         collected_raMOs=zeros(norbs,1);
%         for siteno =   first_site:second_site
%            for orbitalno = first_orbital:last_orbital
%              psi_target = make_target_AO(norbs, siteno, orblistW, orbitalno);
%              fprintf('   RUN %d:  Targets created for site %d orbitals %d:%d.\n',script_run, siteno,first_orbital,last_orbital);
%              [OccupCOs_temp raMOs_H1s DOS_so_far_temp nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%              nraMOs=nraMOs+1;
%              collected_raMOs(:,nraMOs) = raMOs_H1s(:,1);
%            end
%         end
%         movie_on = scriptcodes(script_run,7);   
%         subplot(1,2,2);
%         hold off
%         plot([0 DOSmax],[E_Fermi-Eshift,E_Fermi-Eshift],'color',[0 0 1]);
%         hold on
%         plot_DOS(E_COlist-Eshift,DOS,[.95 .95 .95],[0 0 0],Emin-Eshift, Emax-Eshift, Estep, broadening)
%         xlabel('DOS');
%         ylabel('E [eV]');
%         axis([0 DOSmax Emin-Eshift Emax-Eshift]);
%         hold on
%         axis_prop = gca;
%         axis_prop.Box = 'on';
%         axis_prop.XTick = [0 DOSmax]; 
%          axis_prop.TickDir = 'out';
% 
%         subplot(1,2,1);
%         [DOS_so_far]=raMObands(filename,templatefile,orblistW,H, DOS_scale*DOS_so_far, OccupCOs, COs, collected_raMOs,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%     end
     if (scriptcodes(script_run,1) == 2) 
        first_void =  scriptcodes(script_run,2);
        last_void = scriptcodes(script_run,3);
        rrange= scriptcodes(script_run,4);
         movie_on = scriptcodes(script_run,7);
        for void_no = first_void:last_void
           psi_target = make_target_cluster_sp(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no);
           [psi_previous nelectrons_left]=reconstruct_targetsDFT(filename,templatefile,orblistW, psi_target,movie_on, nelectrons_left, number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states, psi_previous,S_original,templatefile,scalefactor);
        end
     end
%     if (scriptcodes(script_run,1) == 200) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
% %        movie_on = scriptcodes(script_run,7);
%        movie_on = 0;
%        nraMOs=0;
%        DOS_scale = scriptcodes(script_run,6);
%        collected_raMOs=zeros(norbs,1);
%        for void_no = first_void:last_void
%           movie_on = 0;
%           psi_target = make_target_cluster_sp(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no);
%           [OccupCOs_temp raMOs_H1s DOS_so_far_temp nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%           nraMOs=nraMOs+1;
%           collected_raMOs(:,nraMOs) = raMOs_H1s(:,1);
%        end
%        movie_on = scriptcodes(script_run,7);
%        subplot(1,2,2);
%        hold off
%        plot([0 DOSmax],[E_Fermi-Eshift,E_Fermi-Eshift],'color',[0 0 1]);
%        hold on
%        plot_DOS(E_COlist-Eshift,DOS,[.95 .95 .95],[0 0 0],Emin-Eshift, Emax-Eshift, Estep, broadening)
%        xlabel('DOS');
%        ylabel('E [eV]');
%        axis([0 DOSmax Emin-Eshift Emax-Eshift]);
%        hold on
%        axis_prop = gca;
%        axis_prop.Box = 'on';
%        axis_prop.XTick = [0 DOSmax]; 
%  axis_prop.TickDir = 'out';
% 
%        subplot(1,2,1);
%        [DOS_so_far]=raMObands(filename,templatefile,orblistW,H, DOS_scale*DOS_so_far, OccupCOs, COs, collected_raMOs,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%     end    
%     if (scriptcodes(script_run,1) == 22)
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%         movie_on = scriptcodes(script_run,7);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_sp3(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
%     if (scriptcodes(script_run,1) == 23) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%         movie_on = scriptcodes(script_run,7);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_sp3d2(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
%     if (scriptcodes(script_run,1) == 24) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%        first_MO =  scriptcodes(script_run,5);
%        last_MO = scriptcodes(script_run,6);
%        movie_on = scriptcodes(script_run,7);
%        psi_target = make_target_cageMO(H,filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,first_void:last_void);
%        for void_no = first_MO:last_MO
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target(:,void_no),video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
% 
%     if (scriptcodes(script_run,1) == 203) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%        DOS_scale = scriptcodes(script_run,6);
% %        movie_on = scriptcodes(script_run,7);
%        movie_on = 0;
%        nraMOs=0;
%        collected_raMOs=zeros(norbs,1);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_sp3d2(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no);
%           [OccupCOs_temp raMOs_H1s DOS_so_far_temp nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_scale*DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%           nraMOs=nraMOs+1;
%           collected_raMOs(:,nraMOs) = raMOs_H1s(:,1);
%        end
%        movie_on = scriptcodes(script_run,7);
%        subplot(1,2,2);
%        hold off
%        plot([0 DOSmax],[E_Fermi-Eshift,E_Fermi-Eshift],'color',[0 0 1]);
%        hold on
%        plot_DOS(E_COlist-Eshift,DOS,[.95 .95 .95],[0 0 0],Emin-Eshift, Emax-Eshift, Estep, broadening)
%        xlabel('DOS');
%        ylabel('E [eV]');
%        axis([0 DOSmax Emin-Eshift Emax-Eshift]);
%        hold on
%        axis_prop = gca;
%        axis_prop.Box = 'on';
%        axis_prop.XTick = [0 DOSmax]; 
%         axis_prop.TickDir = 'out';
% 
%        subplot(1,2,1);
%        [DOS_so_far]=raMObands(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, collected_raMOs,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%     end
%     if (scriptcodes(script_run,1) == 21)
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%         movie_on = scriptcodes(script_run,7);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_p(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
     if (scriptcodes(script_run,1) == 6) 
        first_void =  scriptcodes(script_run,2);
        last_void = scriptcodes(script_run,3);
        rrange= scriptcodes(script_run,4);
        movie_on = scriptcodes(script_run,7);
        for void_no = first_void:last_void
           psi_target = make_target_hybrid(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no);
           nhybrids=size(psi_target);
           nhybrids=nhybrids(2);
           for hybrid_counter=1:nhybrids
           [psi_previous nelectrons_left]=reconstruct_targetsDFT(filename,templatefile,orblistW, psi_target(:,hybrid_counter),movie_on, nelectrons_left,number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states, psi_previous,S_original,templatefile,scalefactor);
           end
        end
     end
%     if (scriptcodes(script_run,1) == 31) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%        firstMO = scriptcodes(script_run,5);
%        lastMO = scriptcodes(script_run,6);
%        movie_on = scriptcodes(script_run,7);
%        LCcoeff_file = char(scriptfiles(script_run,1));
%        LCcoeff_file = strcat(LCcoeff_file,'-coeff');
%        [LCcoeff] = textread(LCcoeff_file,'%f');
%        LCcoeff = reshape(LCcoeff',lastMO-firstMO+1, lastMO-firstMO+1);
%        for void_no = first_void:last_void
%            psi_target = make_target_cluster_MOs(filename, char(scriptfiles(script_run,1)), norbs, orblistW, rrange,void_no,H);
%            psi_target2 = 0*psi_target;
%            for MOno1=firstMO:lastMO
%                for MOno2=firstMO:lastMO    
%                    psi_target2(:,MOno1) = psi_target2(:,MOno1)+LCcoeff(MOno2,MOno1)*psi_target(:,MOno2);
%                end
%            end
%            for MOno = firstMO:lastMO
%                 [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target2(:,MOno),video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%            end
%        end
%     end
%  %   if (scriptcodes(script_run,1) == 31) 
%  %      first_void =  scriptcodes(script_run,2);
%  %      last_void = scriptcodes(script_run,3);
%  %      rrange= scriptcodes(script_run,4);
%  %      firstMO = scriptcodes(script_run,5);
%  %      lastMO = scriptcodes(script_run,6);
%  %      movie_on = scriptcodes(script_run,7);
%  %      for void_no = first_void:last_void
%  %         fprintf('MOLECULE %d\n',void_no); 
%  %          psi_target = make_target_cluster_MOs_sp_no_pz(filename, char(scriptfiles(script_run,1)), norbs, orblistW, rrange,void_no,firstMO:lastMO,H);
%  %         [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%  %      end
%  %   end
%     if (scriptcodes(script_run,1) == 41)
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%         movie_on = scriptcodes(script_run,7);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_sp2d(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no,1);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
%     if (scriptcodes(script_run,1) == 42) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%         movie_on = scriptcodes(script_run,7);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_sp2d(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no,2);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
%     if (scriptcodes(script_run,1) == 43) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%         movie_on = scriptcodes(script_run,7);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_sp2d_bissect(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no,1);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
%     if (scriptcodes(script_run,1) == 44) 
%        first_void =  scriptcodes(script_run,2);
%        last_void = scriptcodes(script_run,3);
%        rrange= scriptcodes(script_run,4);
%         movie_on = scriptcodes(script_run,7);
%        for void_no = first_void:last_void
%           psi_target = make_target_cluster_sp2d_bissect(filename, char(scriptfiles(script_run,1)), norbs, orblistW,rrange,void_no,2);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left]=reconstruct_targets(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift);
%        end
%     end
%     if (scriptcodes(script_run,1) == 5) 
%           psi_target = 0;
%           movie_on = scriptcodes(script_run,7);
%           [OccupCOs raMOs_H1s DOS_so_far nelectrons_left, psi_previous]=reconstruct_remaindersDFT(filename,templatefile,orblistW,H, DOS_so_far, OccupCOs, COs, psi_target,video_file,movie_on, nelectrons_left, frame_stepsize,E_COlist,DOSparameters,scalefactor,Eshift,number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states, psi_previous,S_original);
%     end
%         
end

end


%%%%%%%%%%%%%%%%%%%%% SUPPORTING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function psi_target = make_target_cluster_sp2d(filename, clusterfile, norbs, orblist, rrange,void_no,mode)
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];


    limit=size(y);
    counter=limit(1);

    % generate sp2d hybrid functions
    h = zeros(9,4);
    h(1,1)=1/2;
    h(1,2)=1/2;
    h(1,3)=1/2;
    h(1,4)=1/2;
    %x-oriented
    h(5,1)=1/2;
    h(2,1)=(1/2)^0.5;
    h(5,2)=1/2;
    h(2,2)=-(1/2)^0.5;
    %y-oriented
    h(5,3)=-1/2;
    h(5,4)=-1/2;
    h(3,3)=(1/2)^0.5;
    h(3,4)=-(1/2)^0.5;
    
    %  h1 -> +x
    %  h2 -> -x
    %  h3 -> +y
    %  h4 -> -y
if(mode==1) 
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
     
    psi_target = zeros(norbs,1);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       [az,el,r_norm]=cart2sph(delta_r(1,1),delta_r(1,2),delta_r(1,3));
                       %based on az decide which hybrid to use 
                       cos_az = cos(az);
                       sin_az = sin(az);
                       h_use = 0;
                       if(abs(cos_az)>abs(sin_az))&&(cos_az > 0)
                           h_use = 1;
                       end
                       if(abs(cos_az)>abs(sin_az))&&(cos_az < 0)
                           h_use = 2;
                       end
                       if(abs(cos_az)<abs(sin_az))&&(sin_az > 0)
                           h_use = 3;
                       end
                       if(abs(cos_az)<abs(sin_az))&&(sin_az < 0)
                           h_use = 4;
                       end
                       if(h_use == 0)
                           fprintf('Improper angle encountered for void %d and atom %d.n',void_no,nx);
                           return;
                       end
                       ao_number = sum(orblist(1:(nx-1)));
                       for j1 = 1:9
                                psi_target(ao_number+j1,1) = h(j1,h_use);
                       end
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);
end
if(mode==2) 
    psi_target = zeros(norbs,5);
    nx = void_no;
    ao_number = sum(orblist(1:(nx-1)));
    psi_target(ao_number+4,1) = 1;
    psi_target(ao_number+6,2) = 1;
    psi_target(ao_number+7,3) = 1;
    psi_target(ao_number+8,4) = 1;
    psi_target(ao_number+9,5) = 1;
end
end

function psi_target = make_target_cluster_sp2d_bissect(filename, clusterfile, norbs, orblist, rrange,void_no,mode)
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];
    limit=size(y);
    counter=limit(1);

    % generate sp2d hybrid functions
    h = zeros(9,4);
    h(1,1)=1/2;
    h(1,2)=1/2;
    h(1,3)=1/2;
    h(1,4)=1/2;
    %++ oriented
    h(7,1)=1/2;
    h(2,1)=(1/2);
    h(3,1)=(1/2); 
    h(7,2)=1/2;
    h(2,2)=-(1/2);
    h(3,2)=-(1/2); 
    %+- oriented
    h(7,3)=-1/2;
    h(7,4)=-1/2;
    h(2,3)=-(1/2);
    h(3,3)=(1/2);
    h(2,4)=(1/2);
    h(3,4)=-(1/2);
    
    %  h1 -> +x +y
    %  h2 -> -x -y
    %  h3 -> -x +y
    %  h4 -> +x -y
    
if(mode==1) 
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');


    psi_target = zeros(norbs,1);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       [az,el,r_norm]=cart2sph(delta_r(1,1),delta_r(1,2),delta_r(1,3));
                       %based on az decide which hybrid to use 
                       cos_az = cos(az);
                       sin_az = sin(az);
                       h_use = 0;
                       if (cos_az > 0)&&(sin_az > 0)
                           h_use = 1;
                       end
                       if (cos_az > 0)&&(sin_az < 0)
                           h_use = 4;
                       end
                       if (cos_az < 0)&&(sin_az > 0)
                           h_use = 3;
                       end
                       if (cos_az < 0)&&(sin_az < 0)
                           h_use = 2;
                       end
                       if(h_use == 0)
                           fprintf('Improper angle encountered for void %d and atom %d.n',void_no,nx);
                           return;
                       end
                       ao_number = sum(orblist(1:(nx-1)));
                       for j1 = 1:9
                                psi_target(ao_number+j1,1) = h(j1,h_use);
                       end
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);
end
if(mode==2) 
    psi_target = zeros(norbs,5);
    nx = void_no;
    ao_number = sum(orblist(1:(nx-1)));
    psi_target(ao_number+4,1) = 1;
    psi_target(ao_number+5,2) = 1;
    psi_target(ao_number+6,3) = 1;
    psi_target(ao_number+8,4) = 1;
    psi_target(ao_number+9,5) = 1;
end
end

function psi_target = make_target_cluster_sp3(filename, clusterfile, norbs, orblist, rrange,void_no)
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];

    limit=size(y);
    counter=limit(1);

    psi_target = zeros(norbs,1);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       ao_number = sum(orblist(1:(nx-1)));
                       psi_target(ao_number+1,1) = (1/2);
                       psi_target(ao_number+2,1) = (3/4)^0.5*delta_r(1,1);
                       psi_target(ao_number+3,1) = (3/4)^0.5*delta_r(1,2);
                       psi_target(ao_number+4,1) = (3/4)^0.5*delta_r(1,3);                       
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);  
end

function psi_target = make_target_cluster_sp3d2(filename, clusterfile, norbs, orblist, rrange,void_no)
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];

    limit=size(y);
    counter=limit(1);

    psi_target = zeros(norbs,1);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       if (delta_r(1,1) > 0) 
                          sp3d2_x(1,1) = 0.5*(2/3)^0.5; 
                          sp3d2_x(2,1) = (1/2)^0.5;
                          sp3d2_x(3,1) = 0.0;
                          sp3d2_x(4,1) = 0.0;
                          sp3d2_x(5,1) = 0.5;
                          sp3d2_x(6,1) = -(1/12)^0.5;
                       else
                          sp3d2_x(1,1) = 0.5*(2/3)^0.5; 
                          sp3d2_x(2,1) = -(1/2)^0.5;
                          sp3d2_x(3,1) = 0.0;
                          sp3d2_x(4,1) = 0.0;
                          sp3d2_x(5,1) = 0.5;
                          sp3d2_x(6,1) = -(1/12)^0.5;                           
                       end
                       if (delta_r(1,2) > 0) 
                          sp3d2_y(1,1) = 0.5*(2/3)^0.5; 
                          sp3d2_y(2,1) = 0;
                          sp3d2_y(3,1) = (1/2)^0.5;
                          sp3d2_y(4,1) = 0.0;
                          sp3d2_y(5,1) = -0.5;
                          sp3d2_y(6,1) = -(1/12)^0.5;
                       else
                          sp3d2_y(1,1) = 0.5*(2/3)^0.5; 
                          sp3d2_y(2,1) = 0.0;
                          sp3d2_y(3,1) = -(1/2)^0.5;
                          sp3d2_y(4,1) = 0.0;
                          sp3d2_y(5,1) = -0.5;
                          sp3d2_y(6,1) = -(1/12)^0.5;                           
                       end
                       if (delta_r(1,3) > 0) 
                          sp3d2_z(1,1) = (1/6)^0.5; 
                          sp3d2_z(2,1) = 0.0;
                          sp3d2_z(3,1) = 0.0;
                          sp3d2_z(4,1) = (1/2)^0.5;
                          sp3d2_z(5,1) = 0.0;
                          sp3d2_z(6,1) = (1/3)^0.5;
                       else
                          sp3d2_z(1,1) = (1/6)^0.5; 
                          sp3d2_z(2,1) = 0.0;
                          sp3d2_z(3,1) = 0.0;
                          sp3d2_z(4,1) = -(1/2)^0.5;
                          sp3d2_z(5,1) = 0.0;
                          sp3d2_z(6,1) = (1/3)^0.5;                           
                       end                      
                       ao_number = sum(orblist(1:(nx-1)));
                       psi_target(ao_number+1,1) = abs(delta_r(1,1))*sp3d2_x(1,1) + abs(delta_r(1,2))*sp3d2_y(1,1)+abs(delta_r(1,3))*sp3d2_z(1,1);
                       psi_target(ao_number+2,1) = abs(delta_r(1,1))*sp3d2_x(2,1) + abs(delta_r(1,2))*sp3d2_y(2,1)+abs(delta_r(1,3))*sp3d2_z(2,1);
                       psi_target(ao_number+3,1) = abs(delta_r(1,1))*sp3d2_x(3,1) + abs(delta_r(1,2))*sp3d2_y(3,1)+abs(delta_r(1,3))*sp3d2_z(3,1);
                       psi_target(ao_number+4,1) = abs(delta_r(1,1))*sp3d2_x(4,1) + abs(delta_r(1,2))*sp3d2_y(4,1)+abs(delta_r(1,3))*sp3d2_z(4,1);
                       psi_target(ao_number+5,1) = abs(delta_r(1,1))*sp3d2_x(5,1) + abs(delta_r(1,2))*sp3d2_y(5,1)+abs(delta_r(1,3))*sp3d2_z(5,1);
                       psi_target(ao_number+6,1) = abs(delta_r(1,1))*sp3d2_x(6,1) + abs(delta_r(1,2))*sp3d2_y(6,1)+abs(delta_r(1,3))*sp3d2_z(6,1);
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);  
end

function psi_target = make_target_cluster_sp(filename, clusterfile, norbs, orblist, rrange,void_no)
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];

    limit=size(y);
    counter=limit(1);

    psi_target = zeros(norbs,1);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       ao_number = sum(orblist(1:(nx-1)));
                       psi_target(ao_number+1,1) = (1/2)^0.5;
                       if(orblist(nx) >1)
                         psi_target(ao_number+2,1) = (1/2)^0.5*delta_r(1,1);
                         psi_target(ao_number+3,1) = (1/2)^0.5*delta_r(1,2);
                         psi_target(ao_number+4,1) = (1/2)^0.5*delta_r(1,3);
                       else
                         fprintf('WARNING:  H 1s orbital being included as cage vertex. \n');
                       end
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);      
end

function psi_target = make_target_cluster_p(filename, clusterfile, norbs, orblist, rrange,void_no)
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];

    limit=size(y);
    counter=limit(1);

    psi_target = zeros(norbs,1);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       ao_number = sum(orblist(1:(nx-1)));
                       psi_target(ao_number+1,1) = 0;
                       psi_target(ao_number+2,1) = delta_r(1,1);
                       psi_target(ao_number+3,1) = delta_r(1,2);
                       psi_target(ao_number+4,1) = delta_r(1,3);                       
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);      
end

function psi_target = make_target_cageMO(H, filename, clusterfile, norbs, orblist, rrange,void_no_range)
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];

    limit=size(y);
    counter=limit(1);
  tcounter = 0;
  nvoids = size(void_no_range);
  nvoids = nvoids(2);
  psi0_target = zeros(norbs,nvoids);
  for void_no = void_no_range
    tcounter = tcounter+1;
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       ao_number = sum(orblist(1:(nx-1)));
                       psi0_target(ao_number+1,tcounter) = (1/2)^0.5;
                       psi0_target(ao_number+2,tcounter) = (1/2)^0.5*delta_r(1,1);
                       psi0_target(ao_number+3,tcounter) = (1/2)^0.5*delta_r(1,2);
                       psi0_target(ao_number+4,tcounter) = (1/2)^0.5*delta_r(1,3);                       
                end
          end
          end
          end
    end
    psi0_target(:,tcounter) = psi0_target(:,tcounter)/norm(psi0_target(:,tcounter));  
  end
  for j = 1:nvoids
      for k = 1:nvoids
      H_mini(j,k) = psi0_target(:,j)'*H*psi0_target(:,k);
      S_mini(j,k) = psi0_target(:,j)'*psi0_target(:,k);
      end
  end
  
  H_mini = (H_mini+H_mini')/2;
  S_mini = (S_mini+S_mini')/2;
  
  [pre_psi,E]=eig(H_mini,S_mini);
  
  psi_target = zeros(norbs,nvoids);
    for j = 1:nvoids
        for k = 1:nvoids
           psi_target(:,j) = psi_target(:,j) + pre_psi(k,j)*psi0_target(:,k);
        end
    end
end

function psi_target = filter_target_cluster_sp(in_psi, filename, clusterfile, norbs, orblist, rrange,void_no)
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];

    limit=size(y);
    counter=limit(1);
    n_raMOs = size(in_psi);
    n_raMOs = n_raMOs(2);
    psi_target = zeros(norbs,n_raMOs);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       ao_number = sum(orblist(1:(nx-1)));
                       for ao_scan = 1:orblist(nx)
                       psi_target(ao_number+ao_scan,raMOno) = in_psi(ao_number+ao_scan,raMOno);
                       end
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);      
end

function psi_target = make_target_cluster_spy(filename, clusterfile, norbs, orblist, rrange,void_no)
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information

    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];

    limit=size(y);
    counter=limit(1);

    psi_target = zeros(norbs,1);
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       ao_number = sum(orblist(1:(nx-1)));
                       psi_target(ao_number+1,1) = (1/2)^0.5;      % s orbital
                       psi_target(ao_number+2,1) = 0;              % px orbital
                       if delta_r(1,2) > 0
                            psi_target(ao_number+3,1) = (1/2)^0.5; % py orbital
                       else
                            psi_target(ao_number+3,1) = -(1/2)^0.5;
                       end    
                       psi_target(ao_number+4,1) = 0;               % pz orbital                
                end
          end
          end
          end
    end
    psi_target = psi_target/norm(psi_target);      
end

function psi_target = make_target_hybrid(filename, clusterfile, norbs, orblist, rrange,void_no)
    [coeff]=textread(clusterfile,'%f');

    ntargets = size(coeff);
    ntargets = ntargets(1)/9;
    psi_target = zeros(norbs,ntargets);
    nx = void_no;
    ao_number = sum(orblist(1:(nx-1)));
    for target_no = 1:ntargets
                       psi_target(ao_number+1,target_no) = coeff((target_no-1)*9+1,1);
                       if(orblist(nx)>1) 
                          psi_target(ao_number+2,target_no) = coeff((target_no-1)*9+2,1);
                          psi_target(ao_number+3,target_no) = coeff((target_no-1)*9+3,1);
                          psi_target(ao_number+4,target_no) = coeff((target_no-1)*9+4,1); 
                       end
                       if(orblist(nx)>4) 
                          psi_target(ao_number+5,target_no) = coeff((target_no-1)*9+5,1);
                          psi_target(ao_number+6,target_no) = coeff((target_no-1)*9+6,1);
                          psi_target(ao_number+7,target_no) = coeff((target_no-1)*9+7,1); 
                          psi_target(ao_number+8,target_no) = coeff((target_no-1)*9+8,1);
                          psi_target(ao_number+9,target_no) = coeff((target_no-1)*9+9,1);
                       end
    end
end

function psi_target = make_target_AO(norbs, atom_list, orblist, orbitals_to_use)
 natoms = size(atom_list);
 natoms = natoms(2);
 norbitals_to_use = size(orbitals_to_use);
 norbitals_to_use = norbitals_to_use(2);
 psi_target = zeros(norbs,1);
 counter=0;
 for siteno = atom_list
    for j = orbitals_to_use
       counter = counter+1;
       ao_number = sum(orblist(1:(siteno-1)))+j;
       psi_target(ao_number,counter) = 1;
    end
 end
end

function psi_target = make_target_massAO(norbs, first_site, second_site, orblist, orbitals_to_use)
 natoms = second_site-first_site+1;
 norbitals_to_use = size(orbitals_to_use);
 norbitals_to_use = norbitals_to_use(2);
 psi_target = zeros(norbs,1);
 counter=0;
 for siteno = first_site:second_site
    for j = orbitals_to_use
       counter = counter+1;
       ao_number = sum(orblist(1:(siteno-1)))+j;
       psi_target(ao_number,counter) = 1;
    end
 end
end

function psi_target = make_target_cluster_MOs(filename, clusterfile, norbs, orblist, rrange,void_no,H)
    fprintf('****  TARGET:  CLUSTER MOs   ****\n');
    [name,void_x, void_y, void_z]=textread(clusterfile,'%s %f %f %f');
    cellfile=strcat(filename,'-cell');
     geofile=strcat(filename,'-geo');
    % load structural information
    nvoids=size(void_x);
    nvoids=nvoids(1);
    [name,x,y,z]=textread(geofile,'%s %f %f %f');
    natoms=size(x);
    natoms=natoms(1);

    a=zeros(3,3);
    [a1,a2,a3]=textread(cellfile,' %f %f %f');
    a=[a1,a2,a3];
    range= [-1,2,-1,2,-1,2];
    fprintf('   Structural information read:  %d centers, %d atoms. \n', nvoids,natoms);

    limit=size(y);
    counter=limit(1);

    modelH = zeros(norbs,norbs);

    %find neighbors coordinating void 
    psi_target = zeros(norbs,1);
    current_cluster_orbital=0;
    for nx = 1:limit(1);
          for j1 = range(1):range(2)
          for k1 = range(3):range(4)
          for l1 = range(5):range(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3); 
                if (xnew-void_x(void_no,1))^2+(ynew-void_y(void_no,1))^2+(znew-void_z(void_no,1))^2 < rrange^2
                       % Neighbor to void found; now add hybrid orbital
                       % pointing to it to the target function
                       delta_r = -[xnew-void_x(void_no,1), ynew-void_y(void_no,1), znew-void_z(void_no,1)];
                       delta_r = delta_r/norm(delta_r);
                       ao_number = sum(orblist(1:(nx-1)));
                       for j2=1:orblist(nx)
                           current_cluster_orbital=current_cluster_orbital+1;
                           psi_target(ao_number+j2,current_cluster_orbital) = 1;
                       end
                end
          end
          end
          end
    end
    fprintf('Number of cluster orbitals:  %d, range = %f \n',current_cluster_orbital,rrange);
    n_cluster_orbitals = current_cluster_orbital;
    H_mini = zeros(n_cluster_orbitals,n_cluster_orbitals);
    for j1=1:n_cluster_orbitals
        for j2=1:n_cluster_orbitals
           H_mini(j1,j2) = psi_target(:,j1)'*H*psi_target(:,j2);  
        end
    end
    H_mini = (H_mini + H_mini')/2;
    [cluster_psi, cluster_E] = eig(H_mini);
    cluster_psi_ao = zeros(norbs,n_cluster_orbitals);
    fprintf('\n');
    for j = 1:n_cluster_orbitals
        fprintf('Cluster MO energy %d:  %f\n',j,cluster_E(j,j));
        for k = 1:n_cluster_orbitals
           cluster_psi_ao(:,j) = cluster_psi_ao(:,j) + cluster_psi(k,j)*psi_target(:,k);
        end
    end
    fprintf('\n');
    psi_target = cluster_psi_ao(:,:);
end     
  

function hamiltonian=readhamiln3(filename,orbnum)
%Draws a band structure from the file filename with 
% containing nbands.  The file should be of the format
% of columns:  band number   energy    occupation
% only the values of the energies matter.
% the function format is bandread(filename,nbands,ymax,yminB)

[j1,j2,j3,orb1,orb2,orb3]=textread(filename,'%s %s %s %f %f %f');
orbcount1=size(orb1);
%orbcouls
nt2=size(orb2)
%orbcount3=size(orb3)
%orbcount4=size(orb4)
%orb4((306-33):306)
%if orbcount1 ~= orbcount2%
%	orb2((orbcount2(1)+1):orbcount1) = zeros((orbcount2(1)+1):orbcount1,1);
%end
%if orbcount1 ~= orbcount3
%        orb3((orbcount3(1)+1):orbcount1) = zeros((orbcount3(1)+1):orbcount1,1);
%end
%if orbcount1 ~= orbcount4
%        orb4((orbcount4(1)+1):orbcount1) = zeros((orbcount4(1)+1):orbcount1,1);
%end%
%
orb1new=reshape(orb1,orbnum,orbcount1(1)/orbnum);
orb2new=reshape(orb2,orbnum,orbcount1(1)/orbnum);
orb3new=reshape(orb3,orbnum,orbcount1(1)/orbnum);
for j=0:(orbcount1(1)/orbnum-1)
	hamiltonian1(1:orbnum,3*j+1)=orb1new(1:orbnum,j+1);
	hamiltonian1(1:orbnum,3*j+2)=orb2new(1:orbnum,j+1);
        hamiltonian1(1:orbnum,3*j+3)=orb3new(1:orbnum,j+1);
end
hamiltonian2=hamiltonian1(1:orbnum,1:orbnum);
hamiltonian=hamiltonian2+transpose(hamiltonian2);
for l=1:orbnum
	hamiltonian(l,l)=.5*hamiltonian(l,l);
end
end


%%%%%%  Vincent Yannello's DFT-raMO Support Functions %%%%%%

function [overlap_target_occupied,E_mat]=calculate_target_occupied_overlap(kpoint_repeating,target_orbitals,number_of_planewaves,G,list_of_kpoints,new_atom_positions,astar,bstar,cstar,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states,number_of_occupied_states,occupied_coefficients,atom_number,N,E,Z,n,spin_up_coefficients,spin_down_coefficients)

switch target_orbitals
case 1
  number_of_orbitals=1;
case 4
  number_of_orbitals=4;
case 9
  number_of_orbitals=9;
end

overlap_target_occupied=zeros(max(number_of_spin_up_states,number_of_spin_down_states),number_of_orbitals,number_of_spin_states);
s_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
p_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
d_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
f_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
px_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
py_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
pz_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
dz2_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
dx2y2_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
dxy_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
dxz_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
dyz_overlap=zeros(number_of_planewaves*number_of_occupied_states,1);
for i=1:number_of_occupied_states
  if kpoint_repeating(i)==0 %If this kpoint has been done before, just copy the values
    if target_orbitals == 1 || target_orbitals==4 || target_orbitals == 9 
      s_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=s_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
    end
    if target_orbitals == 4 || target_orbitals == 9 
      px_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=px_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
      py_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=py_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
      pz_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=pz_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
    end
    if target_orbitals == 9 
      dxy_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=dxy_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
      dyz_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=dyz_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
      dxz_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=dxz_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
      dz2_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=dz2_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
      dx2y2_overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,1)=dx2y2_overlap((i-2)*number_of_planewaves+1:(i-1)*number_of_planewaves,1);
    end
  else %If this kpoint hasn't been done before...
    for j=1:number_of_planewaves
      direction=G(j,:)+list_of_kpoints(i,:);
%      scalingfactor=exp(-2i*pi*(dot(direction,new_atom_positions(atom_number,:))));
      direction=-direction;
      scalingfactor=exp(2i*pi*(dot(direction,new_atom_positions(atom_number,:))));
      direction=direction(1)*astar'+direction(2)*bstar'+direction(3)*cstar';
      %scalingfactor=conj(scalingfactor);
      if target_orbitals == 1 || target_orbitals==4 || target_orbitals == 9 
        if norm(direction)==0;
          s_overlap((i-1)*number_of_planewaves+j,1)=0;
        else
          %Plane wave expansion theroem. Radial portion of wavefunction times the radial component of the expansion
          switch N(1)
          case 1
            s_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*8*Z(1)^2/(norm(direction)^2+Z(1)^2)^2*sqrt(Z(1)*pi);
          case 2
            s_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*8*Z(1)^2*(3*Z(1)^2-norm(direction)^2)/(norm(direction)^2+Z(1)^2)^3*sqrt(Z(1)*pi/3);
          case 3
            s_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*64*Z(1)^4*(Z(1)^2-norm(direction)^2)/(norm(direction)^2+Z(1)^2)^4*sqrt(pi*Z(1)/10);
          case 4
            s_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*96*Z(1)^4*(norm(direction)^4-10*norm(direction)^2*Z(1)^2+5*Z(1)^4)/(norm(direction)^2+Z(1)^2)^5*sqrt(pi*Z(1)/315);
          case 5
            s_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*128/3*Z(1)^6*sqrt(pi*Z(1)/14)*(3*norm(direction)^4-10*norm(direction)^2*Z(1)^2+3*Z(1)^4)/(norm(direction)^2+Z(1)^2)^6;
          case 6
            s_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*128*Z(1)^6*sqrt(pi*Z(1)/462)*(7*Z(1)^6-35*norm(direction)^2*Z(1)^4+21*norm(direction)^4*Z(1)^2-norm(direction)^6)/(norm(direction)^2+Z(1)^2)^7;
          case 7
            s_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*1024*Z(1)^8*sqrt(pi*Z(1)/429)*(Z(1)^6-7*norm(direction)^2*Z(1)^4+7*norm(direction)^4*Z(1)^2-norm(direction)^6)/(norm(direction)^2+Z(1)^2)^8;
          end
        end
      end
      if target_orbitals == 4 || target_orbitals==9 
        if norm(direction)==0;
          p_overlap((i-1)*number_of_planewaves+j,1)=0;
          px_overlap((i-1)*number_of_planewaves+j,1)=0;
          py_overlap((i-1)*number_of_planewaves+j,1)=0;
          pz_overlap((i-1)*number_of_planewaves+j,1)=0;
        else
          switch N(2)
          case 2
            p_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*32i*norm(direction)*Z(2)^3/(norm(direction)^2+Z(2)^2)^3*sqrt(pi*Z(2));
          case 3
            p_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*64i*norm(direction)*Z(2)^3*(5*Z(2)^2-norm(direction)^2)/(norm(direction)^2+Z(2)^2)^4*sqrt(pi*Z(2)/30);
          case 4
            p_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*192i*norm(direction)*Z(2)^5*(5*Z(2)^2-3*norm(direction)^2)/(norm(direction)^2+Z(2)^2)^5*sqrt(pi*Z(2)/105);
          case 5
            p_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*128i/5*norm(direction)*Z(2)^5*(3*norm(direction)^4-42*norm(direction)^2*Z(2)^2+35*Z(2)^4)/(norm(direction)^2+Z(2)^2)^6*sqrt(pi*Z(2)/42);
          case 6
            p_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*1024i/3*norm(direction)*Z(2)^7*(3*norm(direction)^4-14*norm(direction)^2*Z(2)^2+7*Z(2)^4)/(norm(direction)^2+Z(2)^2)^7*sqrt(pi*Z(2)/154);
          end
          %Overlap is calculated for a pz orbital direction along the direction of the planewave, now we need to determine the overlap for the px, py, and pz orbitals in our reference frame.
          if sqrt(direction(1)^2+direction(2)^2)<0.00001 %if the vector lies along the z direction
            cosB=1;
            sinB=0;
            sinA=0;
          else %if the vector doesn't
            cosB=direction(1)/sqrt(direction(1)^2+direction(2)^2);
            sinB=direction(2)/sqrt(direction(1)^2+direction(2)^2);
            sinA=sqrt(direction(1)^2+direction(2)^2)/norm(direction);
          end %if
          cosA=direction(3)/norm(direction);
          %Rotation matrix for p orbitals
          C_1=[sinA*cosB,sinA*sinB,cosA];
          px_overlap((i-1)*number_of_planewaves+j,1)=p_overlap((i-1)*number_of_planewaves+j,1)*C_1(1);
          py_overlap((i-1)*number_of_planewaves+j,1)=p_overlap((i-1)*number_of_planewaves+j,1)*C_1(2);
          pz_overlap((i-1)*number_of_planewaves+j,1)=p_overlap((i-1)*number_of_planewaves+j,1)*C_1(3);
        end
      end
      if target_orbitals == 9 
        if norm(direction)==0;
          d_overlap((i-1)*number_of_planewaves+j,1)=0;
          dx2y2_overlap((i-1)*number_of_planewaves+j,1)=0;
          dz2_overlap((i-1)*number_of_planewaves+j,1)=0;
          dxy_overlap((i-1)*number_of_planewaves+j,1)=0;
          dxz_overlap((i-1)*number_of_planewaves+j,1)=0;
          dyz_overlap((i-1)*number_of_planewaves+j,1)=0;
        else
          temp=4*Z(3)*Z(4)/(Z(3)+Z(4))^2;
          temp=temp^(N(3)+0.5);
          temp=sqrt(n(3)^2+n(4)^2+2*temp*n(3)*n(4));
          temp=1/temp;
          switch N(3)
          case 3
            d_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*(n(3)*temp*(-128*norm(direction)^2*Z(3)^4*sqrt(Z(3)*pi/2)/(norm(direction)^2+Z(3)^2)^4)+n(4)*temp*(-128*norm(direction)^2*Z(4)^4*sqrt(Z(4)*pi/2)/(norm(direction)^2+Z(4)^2)^4));
          case 4
            d_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*(n(3)*temp*(-sqrt(pi*Z(3)/63)*192*Z(3)^4*(7*norm(direction)^2*Z(3)^2-norm(direction)^4)/(norm(direction)^2+Z(3)^2)^5)+n(4)*temp*(-sqrt(pi*Z(4)/63)*192*Z(4)^4*(7*norm(direction)^2*Z(4)^2-norm(direction)^4)/(norm(direction)^2+Z(4)^2)^5));
          case 5
            d_overlap((i-1)*number_of_planewaves+j,1)=scalingfactor*(n(3)*temp*(-sqrt(pi*Z(3)/70)*1024*Z(3)^6*(7*norm(direction)^2*Z(3)^2-3*norm(direction)^4)/(3*(norm(direction)^2+Z(3)^2)^6))+n(4)*temp*(-sqrt(pi*Z(4)/70)*1024*Z(4)^6*(7*norm(direction)^2*Z(4)^2-3*norm(direction)^4)/(3*(norm(direction)^2+Z(4)^2)^6)));
          end
          %Same as for the p orbitals, need to rotate our dz2 orbital onto all of the d orbitals.
          if sqrt(direction(1)^2+direction(2)^2)<0.00001 %if the vector lies along the z direction
            cosB=1;
            sinB=0;
            sinA=0;
          else %if the vector doesn't
            cosB=direction(1)/sqrt(direction(1)^2+direction(2)^2);
            sinB=direction(2)/sqrt(direction(1)^2+direction(2)^2);
            sinA=sqrt(direction(1)^2+direction(2)^2)/norm(direction);
          end %if
          cosA=direction(3)/norm(direction);
          %Useful definitions for the d orbital rotation matrix
          c2B=(cosB*cosB)-(sinB*sinB);
          s2B=2*sinB*cosB;
          %d orbital rotation matrix
          C_2=[sqrt(3)/2*sinA*sinA*c2B,1-3/2*sinA*sinA,sqrt(3)*cosB*sinB*sinA*sinA,sqrt(3)*cosA*sinA*cosB,sqrt(3)*cosA*sinA*sinB];
          dx2y2_overlap((i-1)*number_of_planewaves+j,1)=d_overlap((i-1)*number_of_planewaves+j,1)*C_2(1);
          dz2_overlap((i-1)*number_of_planewaves+j,1)=d_overlap((i-1)*number_of_planewaves+j,1)*C_2(2);
          dxy_overlap((i-1)*number_of_planewaves+j,1)=d_overlap((i-1)*number_of_planewaves+j,1)*C_2(3);
          dxz_overlap((i-1)*number_of_planewaves+j,1)=d_overlap((i-1)*number_of_planewaves+j,1)*C_2(4);
          dyz_overlap((i-1)*number_of_planewaves+j,1)=d_overlap((i-1)*number_of_planewaves+j,1)*C_2(5);
        end
      end
    end
  end
end
%Compiles the energies and overlaps together.
switch target_orbitals
case 1
  overlap=s_overlap;
  E_mat=E(1);
case 4
  overlap=[s_overlap,px_overlap,py_overlap,pz_overlap];
  E_mat=[E(1),E(2),E(2),E(2)];
case 9
  overlap=[s_overlap,px_overlap,py_overlap,pz_overlap,dx2y2_overlap,dz2_overlap,dxy_overlap,dxz_overlap,dyz_overlap];
  E_mat=[E(1),E(2),E(2),E(2),E(3),E(3),E(3),E(3),E(3)];
end
%Creating the actual overlap matrix
if number_of_spin_states==1
  %Creating the actual overlap matrix
  overlap_target_occupied=zeros(number_of_occupied_states,number_of_orbitals);
  for i=1:number_of_occupied_states
      overlap_target_occupied(i,1:number_of_orbitals)=occupied_coefficients(:,i)'*overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,:);
  end
else
  for i=1:number_of_spin_up_states
    overlap_target_occupied(i,:,1)=spin_up_coefficients(:,i)'*overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,:);
  end
  for i=number_of_spin_up_states+1:number_of_spin_up_states+number_of_spin_down_states
    overlap_target_occupied(i-number_of_spin_up_states,:,2)=spin_down_coefficients(:,i-number_of_spin_up_states)'*overlap((i-1)*number_of_planewaves+1:i*number_of_planewaves,:);
  end
end
end

function [vect0,eig0]=sort_eig(matrix0,S0)
%Calculates the eigenvalues and eigenvectors and sorts them based off of the eigenvalues.
[vect0,eig0]=eig(matrix0,S0);
[eig0,ind0]=sort(diag(eig0));
vect0=vect0(:,ind0);
end 

function [type]=outcar_reader
%Reads the OUTCAR file to determine the identity of the atoms
outid=fopen('OUTCAR','r');
found=0;
while ~feof(outid)
  temp=fscanf(outid,'%s', [1]);
  if strcmp(temp,'POTCAR:')
    found=1;
    break;
  end
end
temp=fscanf(outid,'%s', [1]);
temp=fscanf(outid,'%s', [1]);
type{found}=temp;
temp=fscanf(outid,'%s', [1]);
stopage=0;
while stopage==0
  temp=fscanf(outid,'%s', [1]);
  temp=fscanf(outid,'%s', [1]);
  temp=fscanf(outid,'%s', [1]);
  for i=1:found
    if strcmp(type{i},temp)
      stopage=1;
    end
  end
  if stopage==1
     break;
  end
  found=found+1;
  type{found}=temp;
  if found==5
    fprintf('ERROR reading OUTCAR for atom identities!\n');
    break;
  end
  temp=fscanf(outid,'%s', [1]);
end
fclose(outid);
end

function atomnumbers=atomlookup(atomidentities)
%Assigns atomic numbers based off of the atomic symbol
numatoms=size(atomidentities,2);
atomnumbers=zeros(1,numatoms);
for i=1:numatoms
  switch atomidentities{i}
  case 'H'
    atomnumbers(i)=1;
  case 'He'
    atomnumbers(i)=2;
  case 'Li'
    atomnumbers(i)=3;
  case 'Be'
    atomnumbers(i)=4;
  case 'B'
    atomnumbers(i)=5;
  case 'C'
    atomnumbers(i)=6;
  case 'N'
    atomnumbers(i)=7;
  case 'O'
    atomnumbers(i)=8;
  case 'F'
    atomnumbers(i)=9;
  case 'Ne'
    atomnumbers(i)=10;
  case 'Na'
    atomnumbers(i)=11;
  case 'Mg'
    atomnumbers(i)=12;
  case 'Al'
    atomnumbers(i)=13;
  case 'Si'
    atomnumbers(i)=14;
  case 'P'
    atomnumbers(i)=15;
  case 'S'
    atomnumbers(i)=16;
  case 'Cl'
    atomnumbers(i)=17;
  case 'Ar'
    atomnumbers(i)=18;
  case 'K'
    atomnumbers(i)=19;
  case 'Ca'
    atomnumbers(i)=20;
  case 'Sc'
    atomnumbers(i)=21;
  case 'Ti'
    atomnumbers(i)=22;
  case 'V'
    atomnumbers(i)=23;
  case 'Cr'
    atomnumbers(i)=24;
  case 'Mn'
    atomnumbers(i)=25;
  case 'Fe'
    atomnumbers(i)=26;
  case 'Co'
    atomnumbers(i)=27;
  case 'Ni'
    atomnumbers(i)=28;
  case 'Cu'
    atomnumbers(i)=29;
  case 'Zn'
    atomnumbers(i)=30;
  case 'Ga'
    atomnumbers(i)=31;
  case 'Ge'
    atomnumbers(i)=32;
  case 'As'
    atomnumbers(i)=33;
  case 'Se'
    atomnumbers(i)=34;
  case 'Br'
    atomnumbers(i)=35;
  case 'Kr'
    atomnumbers(i)=36;
  case 'Rb'
    atomnumbers(i)=37;
  case 'Sr'
    atomnumbers(i)=38;
  case 'Y'
    atomnumbers(i)=39;
  case 'Zr'
    atomnumbers(i)=40;
  case 'Nb'
    atomnumbers(i)=41;
  case 'Mo'
    atomnumbers(i)=42;
  case 'Tc'
    atomnumbers(i)=43;
  case 'Ru'
    atomnumbers(i)=44;
  case 'Rh'
    atomnumbers(i)=45;
  case 'Pd'
    atomnumbers(i)=46;
  case 'Ag'
    atomnumbers(i)=47;
  case 'Cd'
    atomnumbers(i)=48;
  case 'In'
    atomnumbers(i)=49;
  case 'Sn'
    atomnumbers(i)=50;
  case 'Sb'
    atomnumbers(i)=51;
  case 'Te'
    atomnumbers(i)=52;
  case 'I'
    atomnumbers(i)=53;
  case 'Xe'
    atomnumbers(i)=54;
  case 'Cs'
    atomnumbers(i)=55;
  case 'Ba'
    atomnumbers(i)=56;
  case 'La'
    atomnumbers(i)=57;
  case 'Ce'
    atomnumbers(i)=58;
  case 'Pr'
    atomnumbers(i)=59;
  case 'Nd'
    atomnumbers(i)=60;
  case 'Pm'
    atomnumbers(i)=61;
  case 'Sm'
    atomnumbers(i)=62;
  case 'Eu'
    atomnumbers(i)=63;
  case 'Gd'
    atomnumbers(i)=64;
  case 'Tb'
    atomnumbers(i)=65;
  case 'Dy'
    atomnumbers(i)=66;
  case 'Ho'
    atomnumbers(i)=67;
  case 'Er'
    atomnumbers(i)=68;
  case 'Tm'
    atomnumbers(i)=69;
  case 'Yb'
    atomnumbers(i)=70;
  case 'Lu'
    atomnumbers(i)=71;
  case 'Hf'
    atomnumbers(i)=72;
  case 'Ta'
    atomnumbers(i)=73;
  case 'W'
    atomnumbers(i)=74;
  case 'Re'
    atomnumbers(i)=75;
  case 'Os'
    atomnumbers(i)=76;
  case 'Ir'
    atomnumbers(i)=77;
  case 'Pt'
    atomnumbers(i)=78;
  case 'Au'
    atomnumbers(i)=79;
  case 'Hg'
    atomnumbers(i)=80;
  case 'Tl'
    atomnumbers(i)=81;
  case 'Pb'
    atomnumbers(i)=82;
  case 'Bi'
    atomnumbers(i)=83;
  case 'Po'
    atomnumbers(i)=84;
  case 'At'
    atomnumbers(i)=85;
  case 'Rn'
    atomnumbers(i)=86;
  case 'Fr'
    atomnumbers(i)=87;
  case 'Ra'
    atomnumbers(i)=88;
  case 'Ac'
    atomnumbers(i)=89;
  case 'Th'
    atomnumbers(i)=90;
  case 'Pa'
    atomnumbers(i)=91;
  case 'U'
    atomnumbers(i)=92;
  case 'Np'
    atomnumbers(i)=93;
  case 'Pu'
    atomnumbers(i)=94;
  case 'Am'
    atomnumbers(i)=95;
  case 'Cm'
    atomnumbers(i)=96;
  case 'Bk'
    atomnumbers(i)=97;
  case 'Cf'
    atomnumbers(i)=98;
  case 'Es'
    atomnumbers(i)=99;
  case 'Fm'
    atomnumbers(i)=100;
  case 'Md'
    atomnumbers(i)=101;
  case 'No'
    atomnumbers(i)=102;
  case 'Lr'
    atomnumbers(i)=103;
  case 'Rf'
    atomnumbers(i)=104;
  case 'Db'
    atomnumbers(i)=105;
  case 'Sg'
    atomnumbers(i)=106;
  case 'Bh'
    atomnumbers(i)=107;
  case 'Hs'
    atomnumbers(i)=108;
  case 'Mt'
    atomnumbers(i)=109;
  case 'Ds'
    atomnumbers(i)=110;
  case 'Rg'
    atomnumbers(i)=111;
  case 'Cn'
    atomnumbers(i)=112;
  case 'Nh'
    atomnumbers(i)=113;
  case 'Fl'
    atomnumbers(i)=114;
  case 'Mc'
    atomnumbers(i)=115;
  case 'Lv'
    atomnumbers(i)=116;
  case 'Ts'
    atomnumbers(i)=117;
  case 'Og'
    atomnumbers(i)=118;
  end
end
end

function [numatoms,symbol,atom_type,number_of_electrons,number_of_orbitals,N,E,Z1,Z2,n1,n2]=readeht_parms
%Gets the Huckel Parameters from the file below
fileID = fopen('DFT_raMO_eht_parms.dat','r');
strlist=[];
while ~feof(fileID)
   temp=fgetl(fileID);
   if ~strncmp(temp,';',1)&&~isempty(strtrim(temp)) strlist=[strlist,temp,char(10)]; end
end
temp=textscan(strlist,'%s %f %f %f %f %s %f %f %f %f %f %*[^\n]');
atomlist=[1,find(diff(temp{2}'))+1];
numatoms=numel(atomlist);
symbol=temp{1}(atomlist)';
number_of_electrons=temp{3}(atomlist)';
atom_type=temp{2}(atomlist)';
number_of_orbitals=diff([atomlist,numel(temp{2})+1]);
N=zeros(4,numatoms);
E=zeros(4,numatoms);
Z1=zeros(4,numatoms);
Z2=zeros(4,numatoms);
n1=zeros(4,numatoms);
n2=zeros(4,numatoms);
for j=1:numatoms
   N(1:number_of_orbitals(j),j)=temp{5}(atomlist(j):atomlist(j)+number_of_orbitals(j)-1);
   E(1:number_of_orbitals(j),j)=temp{7}(atomlist(j):atomlist(j)+number_of_orbitals(j)-1);
   Z1(1:number_of_orbitals(j),j)=temp{8}(atomlist(j):atomlist(j)+number_of_orbitals(j)-1);
   Z2(1:number_of_orbitals(j),j)=temp{9}(atomlist(j):atomlist(j)+number_of_orbitals(j)-1);
   n1(1:number_of_orbitals(j),j)=temp{10}(atomlist(j):atomlist(j)+number_of_orbitals(j)-1);
   n2(1:number_of_orbitals(j),j)=temp{11}(atomlist(j):atomlist(j)+number_of_orbitals(j)-1);
end
number_of_orbitals=number_of_orbitals.^2;
fclose(fileID);
end

function [number_of_electrons,number_of_orbitals,N,E,Z,n]=geteht_parms(atomsymbol,eht_symbol, eht_number_of_orbitals, eht_number_of_electrons, eht_N, eht_E, eht_Z1, eht_Z2, eht_n1, eht_n2)
%Loop until we find the correct atom
i=find(strcmpi(atomsymbol,eht_symbol));

if isempty(i)
   error(['Could not find element ',atomsymbol{1},' in eht_parms.dat file.']);
else %if found, set parameters
   number_of_orbitals=eht_number_of_orbitals(i);
   number_of_electrons=eht_number_of_electrons(i);
   N=eht_N(:,i);
   E=eht_E(:,i);
   Z(1)=eht_Z1(1,i);
   Z(2)=eht_Z1(2,i);
   Z(3)=eht_Z1(3,i);
   Z(4)=eht_Z2(3,i);
   Z(5)=eht_Z1(4,i);
   Z(6)=eht_Z2(4,i);
   n(1)=eht_n1(1,i);
   n(2)=eht_n1(2,i);
   n(3)=eht_n1(3,i);
   n(4)=eht_n2(3,i);
   n(5)=eht_n1(4,i);
   n(6)=eht_n2(4,i);
end
end


function S=make_overlap_mat(occupied_coefficients,kpoint_repeating)
%Makes the overlap matrix from the occupied states.
number_of_occupied_states=size(occupied_coefficients,2);
S=zeros(number_of_occupied_states,number_of_occupied_states);
for i=1:size(kpoint_repeating,1)
  if kpoint_repeating(i)==1 && i>1
    S(ranger,ranger)=occupied_coefficients(:,ranger)'*occupied_coefficients(:,ranger);
    ranger=[i];
  elseif i==1
    ranger=[i];
  elseif i==size(kpoint_repeating,1)
    ranger=[ranger,i];
    S(ranger,ranger)=occupied_coefficients(:,ranger)'*occupied_coefficients(:,ranger);
  else
    ranger=[ranger,i];
  end
end
end

function kpoint=read_kpoints
%Reads in the Kpoints from the KPOINTS file.
fin = fopen('KPOINTS','r');
input_name=fscanf(fin,'%s',[1]);
centering=fscanf(fin,'%i',[1]);
centername=fscanf(fin,'%s',[1]);
kpoint=fscanf(fin,'%i %i %i',[3]);
shift=fscanf(fin,'%f %f %f',[3]);
end


function make_readableparams(input_name,name_of_readableparams,occupied_filling_level,use_fermi_energy,Emin,Emax)
%Takes the output from WaveTrans90 and makes it easier to read by MATLAB.
fout = fopen(name_of_readableparams,'w');
fileID = fopen(input_name,'r');
%Reads in the Wavefunction Information from the Input File
number_of_spin_states=fscanf(fileID,'%i',1);
number_of_kpoints=fscanf(fileID,'%i',1);
number_of_bands=fscanf(fileID,'%i',1);
a=fscanf(fileID,'%f %f %f',3);
b=fscanf(fileID,'%f %f %f',3);
c=fscanf(fileID,'%f %f %f',3);
astar=fscanf(fileID,'%f %f %f',3);
bstar=fscanf(fileID,'%f %f %f',3);
cstar=fscanf(fileID,'%f %f %f',3);
list_of_kpoints=zeros(1,3);

%list_of_kpoints(1,:)=fscanf(fileID,'%f %f %f', 3);
%current_band=fscanf(fileID,'%i',1);
%number_of_planewaves=fscanf(fileID,'%i',1);
%G=zeros(number_of_planewaves,3);
%occupied_coefficients=zeros(number_of_planewaves,1);
%temp=fgetl(fileID);
%temp=fgetl(fileID);
%strlist=[];
%strlist=[strlist,temp,char(10)];
%temp=textscan(strlist,'%c %f %c %f %c %f');
%Energies=temp{2};
%Filling=temp{6};
%for j=1:number_of_planewaves
%  temp=fgetl(fileID);
%  strlist=[];
%  strlist=[strlist,temp,char(12)];
%  temp=textscan(strlist,'%f %f %f %c %f %c %f %c');
%  G(j,:)=[temp{1},temp{2},temp{3}];
%  occupied_coefficients(j,1)=complex(temp{5},temp{7});
%end
G = [0 0 0];
number_of_planewaves=1;
if number_of_spin_states==2
  number_of_spin_up_states=0;
  number_of_spin_down_states=0;
end
number_of_occupied_states=0;
for m=1:number_of_spin_states
  fprintf('Reading spin state %d...\n',m);
  for l=1:number_of_kpoints %Looping through all of the kpoints
    foundit=0;
    ktemp=fscanf(fileID,'%f %f %f',3);
    occupied_states_by_k = 0;
    for i=1:number_of_bands
        current_band=fscanf(fileID,'%i',1);
        numplanewavestemp=fscanf(fileID,'%i',1);
        temp=fgetl(fileID);
        temp=fgetl(fileID);
        strlist=[];
        strlist=[strlist,temp,char(20)];
        temp=textscan(strlist,'%c %f %c %f %c %f');
        if (size(temp{2})>1)==[0, 0]
          Energies=temp{2};
          Filling=temp{6};
        else
          Energies=temp{2}(1);
          Filling=temp{6}*10^(temp{2}(2));
        end
        %Will only read the information for occupied states      
        if (Energies >= Emin && Energies <= Emax)
          if number_of_spin_states==2
            if m==1
              number_of_spin_up_states=number_of_spin_up_states+1;
            else
              number_of_spin_down_states=number_of_spin_down_states+1;
            end
          end
          occupied_states_by_k = occupied_states_by_k +1;
          number_of_occupied_states=number_of_occupied_states+1;
          k=1;
          list_of_kpoints(number_of_occupied_states,:)=ktemp;
          if(number_of_occupied_states)==1
             kpoint_repeating=1;
          else
            if list_of_kpoints(number_of_occupied_states,:)==list_of_kpoints(number_of_occupied_states-1,:)
              temp=[kpoint_repeating;0];
            else
              temp=[kpoint_repeating;1];
            end
            kpoint_repeating=temp;
          end
          for j=1:numplanewavestemp
            temp=fgetl(fileID);
            strlist=[];
            strlist=[strlist,temp,char(12)];
            temp=textscan(strlist,'%f %f %f %c %f %c %f %c');
            found=0;
            start=k;
 %           if(number_of_planewaves > 0)
              while found==0
                if G(k,:)==[temp{1},temp{2},temp{3}]
                  found=1;
                  break;
                end
                if k==number_of_planewaves
                  k=0;
                end
                k=k+1;
                if k==start
                  break;
                end
              end
              if found==1
                occupied_coefficients(k,number_of_occupied_states)=complex(temp{5},temp{7});
              else
                number_of_planewaves=number_of_planewaves+1;
                G(number_of_planewaves,:)=[temp{1},temp{2},temp{3}];
                occupied_coefficients(number_of_planewaves,number_of_occupied_states)=complex(temp{5},temp{7});
              end
%            else
%                number_of_planewaves=number_of_planewaves+1;
%                G(number_of_planewaves,:)=[temp{1},temp{2},temp{3}];
%                occupied_coefficients(number_of_planewaves,number_of_occupied_states)=complex(temp{5},temp{7});
%            end
          end
        else
          if foundit==0
 %           fprintf('For K-Point Number %i, Spin State %i, found %i Occupied States\n', l, m, i-1);
            foundit=1;
          end
          for j=1:numplanewavestemp
            temp=fgetl(fileID);
          end
        end
    end
    fprintf('For K-Point Number %i, Spin State %i, found %i Occupied States\n', l, m, occupied_states_by_k);
  end
end
if number_of_spin_states==1
  number_of_spin_up_states=0;
  number_of_spin_down_states=0;
end
fprintf(fout,'%i\n', number_of_spin_states);
fprintf(fout,'%i\n', number_of_spin_up_states);
fprintf(fout,'%i\n', number_of_spin_down_states);
fprintf(fout,'%i\n', number_of_occupied_states);
fprintf(fout,'%i\n', number_of_planewaves);
fprintf(fout,'%f %f %f\n', list_of_kpoints');
fprintf(fout,'%i %i %i\n', G');
for i=1:number_of_occupied_states
  temp=[real(occupied_coefficients(:,i)), imag(occupied_coefficients(:,i))];
  fprintf(fout,'%.20f %.20f\n', temp');
end
fprintf('Made readableparams file!\n');
end

function   [number_of_occupied_states,number_of_planewaves,list_of_kpoints,G,occupied_coefficients,kpoint_repeating,number_of_spin_states,number_of_spin_up_states,number_of_spin_down_states]=read_readableparams(name_of_readableparams)
%Reads in the readableparams.txt file.
fin = fopen(name_of_readableparams,'r');
number_of_spin_states=fscanf(fin,'%i',1);
number_of_spin_up_states=fscanf(fin,'%i',1);
number_of_spin_down_states=fscanf(fin,'%i',1);
number_of_occupied_states=fscanf(fin,'%i',1);
number_of_planewaves=fscanf(fin,'%i',1);
kpoint_repeating=zeros(number_of_occupied_states,1);
kpoint_repeating(1)=1;
list_of_kpoints=fscanf(fin,'%f %f %f',[3,number_of_occupied_states])';
for i=2:number_of_occupied_states
  if list_of_kpoints(i,:)==list_of_kpoints(i-1,:)
  else
    kpoint_repeating(i)=1;
  end
end
G=fscanf(fin,'%i %i %i',[3,number_of_planewaves])';
occupied_coefficients=zeros(number_of_planewaves,number_of_occupied_states);
for i=1:number_of_occupied_states
  for j=1:number_of_planewaves
    temp1=fscanf(fin,'%f %f',[2,1]);
    occupied_coefficients(j,i)=complex(temp1(1),temp1(2));
  end
end
fprintf('Read readableparams file!\n');
end

function [acell,bcell,ccell,number_of_atoms_by_type,total_number_of_atoms,atom_positions]=read_poscar(number_of_atom_types)
%Reads in the POSCAR information.

f3=fopen('POSCAR','r');
temp=fscanf(f3,'%s',1);
fprintf('%s \n',temp);
temp=fscanf(f3,'%f',1);
acell=fscanf(f3,'%f %f %f',3);
bcell=fscanf(f3,'%f %f %f',3);
ccell=fscanf(f3,'%f %f %f',3);
number_of_atoms_by_type=zeros(number_of_atom_types,1);
for i=1:number_of_atom_types
  number_of_atoms_by_type(i,1)=fscanf(f3,'%i',1);
end
total_number_of_atoms=sum(number_of_atoms_by_type);
temp=fscanf(f3,'%s',1);
atom_positions=zeros(total_number_of_atoms,3);
fprintf('# of atoms in POSCAR: %i \n', total_number_of_atoms);
for i=1:total_number_of_atoms
  atom_positions(i,1)=fscanf(f3,'%f',1);
  atom_positions(i,2)=fscanf(f3,'%f',1);
  atom_positions(i,3)=fscanf(f3,'%f',1);
end
fclose(f3);
end


