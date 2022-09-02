%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                          UNIVARIATE REGRESSION MODEL                         %
%                                                                              %
%                                 January 2019                                 %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program has been tested with Matlab R2015a

% ==============================================================================
% Instructions:

% Edit parameters marked with an arrow like this: <---
 
% ==============================================================================
% POLYNOMIAL FIT:

             x = []; % <--- independent variable
             y = []; % <--- dependent variable

centerscale_sw = 1;  % <--- center and rescale (1 = yes, 0 = no)
  safecoeff_sw = 1;  % <--- 1 = compute only coefficients that produce no warnings, 0 = go bananas

clc
clearvars p S n delta all_delta all_delta_sum all_predicted_y
for n = 1:length(x);
    lastwarn('')
    if centerscale_sw
        [p,S,mu] = polyfit(x,y,n);
    else
        [p,S] = polyfit(x,y,n);
    end
    if safecoeff_sw && ~isempty(lastwarn)
        delta = Inf;
        this_predicted_y = NaN;
        clc
    else
        all_p{n} = p;
        if centerscale_sw
            [this_predicted_y,delta] = polyval(p,x,S,mu);
        else
            [this_predicted_y,delta] = polyval(p,x,S);
        end
    end
    clc
    all_predicted_y(n,:) = this_predicted_y;
    all_delta(n,:) = delta;
    all_delta_sum(n) = sum(delta(:));
end
best_n = find(all_delta_sum == min(all_delta_sum));

disp('all predicted y = ')
disp(all_predicted_y)

% %% .............................................................................
% DOUBLE-CHECK AND MAKE EQUATION:

chosen_n = best_n; % <--- choose polynomial order to build equation
% chosen_n = 3;

% double-check:
clearvars check_predicted_y check_delta
check_predicted_y = zeros(1,length(x));
check_delta = check_predicted_y;
for i_0 = 1:length(x)
    check_predicted_y(i_0) = 0;
    for i = 1:chosen_n+1
        check_predicted_y(i_0) = check_predicted_y(i_0) + all_p{chosen_n}(i) * x(i_0)^(chosen_n-(i-1));
    end
    check_delta(i_0) = all_predicted_y(chosen_n,i_0) - check_predicted_y(i_0);
end

% make equation;
equation_str = '';
sign_str = '';
chosen_p = all_p{chosen_n};
abs_chosen_p = abs(chosen_p);
for i = 1:chosen_n+1
    if chosen_n-(i-1) == 0
        power_str = ')';
    elseif chosen_n-(i-1) == 1
        power_str = ' * x)';
    else
        power_str = sprintf(' * x^%g)',chosen_n-(i-1));
    end
    if i > 1
        if sign(chosen_p(i)) == 1
            sign_str = '+';
        elseif sign(chosen_p(i)) == -1
            sign_str = '-';
        end
    end
    equation_str = sprintf('%s%s (%.3g%s ',equation_str,sign_str,abs_chosen_p(i),power_str);
end
disp(sprintf('  best polynomial order n = %i',best_n))
disp(' ')
disp(sprintf('chosen polynomial order n = %i',chosen_n))
disp(' ')
disp('equation: ')
disp(sprintf('y =%s',equation_str))
disp(' ')
disp(sprintf('polyfit predicted at chosen n, y = %s',mat2str(all_predicted_y(chosen_n,:),3)))
disp(sprintf('   loop predicted at chosen n, y = %s',mat2str(check_predicted_y,3)))
disp(sprintf('                      difference = %s',mat2str(check_delta,3)))
