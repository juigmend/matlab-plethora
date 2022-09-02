%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                   MAGNITUDE                                  %
%                                                                              %
%                                 October, 2020                                %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyväskylä                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tested on Octave 4

%===============================================================================

% DESCRIPTION:

%   Test different methods for computing magnitude of data.

% INSTRUCTIONS:

%   1) Modify the values indicated by an arrow (<---)
%   2) Run the script and enjoy.

%===============================================================================

clc
clear
close all
tic

%===============================================================================

n_cases = 60;   % <--- number of cases
  n_obs = 30;   % <--- number of observations

    n_d = 4;    % <--- differences
  
   n_MC = 1;  % <--- number of montecarlo realisations
 i_MC_p = 1;    % <--- index of realisations to plot

%-------------------------------------------------------------------------------

for i_MC = 1:n_MC

  X = cell(n_cases,1);

  for i_case = 1:n_cases

    X{i_case} = rand(n_obs); % random data
    
  end  

  X_d = cell(n_d,1);

  % differences (derivatives):
  for i_d = 1:n_d
   
    X_d{i_d} = cell(n_cases,1);
   
    for i_case = 1:n_cases
    
      X_d{i_d}{i_case} = diff( X{i_case} , i_d ); % difference
    end
  end
      
      
  % features from derivatives:
  for i_d = 1:n_d
   
    for i_case = 1:n_cases

      dfeats.L1norm_dpw1{i_d}(i_case) = norm( X_d{i_d}{i_case} , 1 ); %  L1 norm (Manhattan distance)
    
      dfeats.L2norm{i_d}(i_case) = norm( X_d{i_d}{i_case} , 2 ); % L2 norm (Euclidean distance)
      
      dfeats.abssum_B{i_d}(i_case) = sum( abs( X_d{i_d}{i_case}(:)  )); % sum of absolute values
    end
    
  end

  lbl_dfeats = fieldnames(dfeats);
  n_dfeats = length(lbl_dfeats);

  n_sp = n_dfeats * n_d;
    
  for i_dfeats = 1:n_dfeats

    for i_d = 1:n_d
    
      M{i_d}(:,i_dfeats) = dfeats.( lbl_dfeats{i_dfeats}  ){i_d}; 
    
    end
  end


  for i_d = 1:n_d
    
    r{i_d} = corr(M{i_d});
    
    if i_MC == 1
    
      sum_r{i_d} = r{i_d};
    else
    
      sum_r{i_d} = sum_r{i_d} + r{i_d};
    end
  end
  
  
  if any( i_MC == i_MC_p )
    
    i_sp = 0;
      
    for i_dfeats = 1:n_dfeats

      for i_d = 1:n_d
      
        i_sp = i_sp + 1;
        subplot(n_dfeats,n_d,i_sp)
        plot( M{i_d}(:,i_dfeats) )
        
        this_featlbl = lbl_dfeats{i_dfeats};
        this_featlbl( this_featlbl == '_') = ' ';
        title( sprintf( '%s d = %i ', this_featlbl , i_d ) )
      end
    end

    figure

    for i_d = 1:n_d

      lbl_d = num2str(i_d);

      this_corr_lbl = sprintf( 'abs corr d = %i', i_d);
      
      subplot(1,n_d,i_d)
      imagesc(abs( r{i_d} ))
      title(this_corr_lbl)
      
      disp(' ')
      disp(this_corr_lbl)
      disp(r{i_d})
    end
  end
end

if n_MC > 1
    
    figure
    
    for i_d = 1:n_d
        
        lbl_d = num2str(i_d);
        
        this_corr_lbl = sprintf( 'MC sum abs corr d = %i', i_d);
        
        subplot(1,n_d,i_d)
        imagesc( sum_r{i_d} )
        title(this_corr_lbl)
    end
end

toc