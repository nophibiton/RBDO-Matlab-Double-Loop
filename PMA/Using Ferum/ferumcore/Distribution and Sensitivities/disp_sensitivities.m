function d = disp_sensitivities (lsf,formresults,probdata,gfundata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routines display the sensitivity analysis results corresponding         %
% to the analysis of a given limit-state function                              %
%                                                                              %
% In order that the information given by this subroutine are correct the       %
% call to this subroutine must follow a the form analysis of the limit-state   %
% function you want to analyse. This particularly important when performing    %
% a system analysis involving several limit-state function                     %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display SENSITIVITIES WITH RESPECT TO THE DISTRIBUTION PARAMETERS

nrv = size(probdata.marg,1);
beta_sensi = formresults.beta_sensi_thetaf;
pf_sensi = formresults.pf_sensi_thetaf;

   
   disp(sprintf('SENSITIVITIES OF THE RELIABILITY INDEX WITH RESPECT TO DISTRIBUTION PARAMETERS'))
   disp(sprintf('----------------------------------------------------------------------------------------------'))
   disp(sprintf(' var          mean        std dev           par1           par2           par3           par4'))
   for k = 1:nrv
    disp(sprintf('%3.d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n',k,beta_sensi(k,1),beta_sensi(k,2),beta_sensi(k,3),beta_sensi(k,4),beta_sensi(k,5),beta_sensi(k,6)))
   end
   disp(sprintf('----------------------------------------------------------------------------------------------\n'))
   
   disp(sprintf('SENSITIVITIES of THE FAILURE PROBABILITY WITH RESPECT TO DISTRIBUTION PARAMETERS'))
   disp(sprintf('----------------------------------------------------------------------------------------------'))
   disp(sprintf(' var          mean        std dev           par1           par2           par3           par4'))
   for k = 1:nrv
    disp(sprintf('%3.d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n',k,pf_sensi(k,1),pf_sensi(k,2),pf_sensi(k,3),pf_sensi(k,4),pf_sensi(k,5),pf_sensi(k,6)))
   end
   disp(sprintf('----------------------------------------------------------------------------------------------\n'))
   
% Display SENSITIVITIES WITH RESPECT TO THE LIMIT-STATE FUNCTION PARAMETERS
switch lower(gfundata(lsf).parameter)
   case 'yes' 
       switch lower(gfundata(lsf).evaluator)
           case 'basic'
               n_thetag = length(gfundata(lsf).dgthetag);
           otherwise
               n_thetag = 1;
       end
       beta_sensi_thetag = formresults.beta_sensi_thetag;
       pf_sensi_thetag =  formresults.pf_sensi_thetag;
       disp(sprintf('SENSITIVITIES WITH RESPECT TO THE LIMIT-STATE FUNCTION PARAMETERS'))
       disp(sprintf('----------------------------------------------------------------------------------------------'))
       disp(sprintf(' par   d(beta)/d(par)    d(pf1)/d(par)  '))
       for k = 1:n_thetag
        disp(sprintf('%3.d%17.5e%17.5e',k,beta_sensi_thetag(k),pf_sensi_thetag(k)))
       end
       disp(sprintf('----------------------------------------------------------------------------------------------\n'))   
   case 'no'
   otherwise
   disp('Unknown type of limit-state function parameter input.')
end
   
   