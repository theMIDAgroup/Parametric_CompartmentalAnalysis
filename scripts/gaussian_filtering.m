function [B]=gaussian_filtering(A,filter_size,sigma)
% [B]=gaussian_filtering(A,filter_size,sigma)) computes the filtered matrix B 
% with a Gaussian filter to remove Poisson noise on the original image A.

% A is the original image to filter
% filter_size is the dimension of the gaussian window (it has to be an odd number)
% sigma is the standard deviation of the gaussian filter

B=zeros(size(A)); % inizialization matrix B post gaussian filtering
L=filter_size;

doit=false;
if mod(L,2)~=0 % control L=filter_size odd number
    doit=true;
else
    disp('WARNING function ''gaussian_filtering''!')
    disp('incorrect input: [B]=gaussian_filtering(A,filter_size,sigma)')
    disp('filter_size has to be an odd number')
end

% Compute the 2D Gaussian filter with standard deviation sigma (and mean zero)
N=(L-1)/2;
[x, y]=meshgrid(-N:N,-N:N);
g=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
% Normalize to sum to one 
g=g./sum(g(:)); 

% otherwise
% g = 1/(2*pi*sigma^2)*exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));

if doit==true % if the control on the window filter size is successful
    
    n_plus=(L-1)/2; % rows/columns to add on the original matrix
    A_rep=repmat(A,3,3); % replicate the matrix A in all the directions
    % <--> 3 times for rows, 3 times for columns
    
    % To isolate the central 'A' matrix plus (L-1)/2 rows and (L-1)/2 
    % columns (added side by side), compute the number of rows/columns 
    % to eliminate from the replicated matrix 'A_rep'
    waste_rows=size(A,1)-n_plus+1;
    waste_columns=size(A,2)-n_plus+1;
    
    % Select the part of the replicated matrix 'A_rep' on which the filter should acts
    A_big=A_rep(waste_rows:size(A_rep,1)-waste_rows+1, waste_columns:size(A_rep,2)-waste_columns+1);
    
    for i=1:size(A_big,1)-L+1     %| ==> moove the window on all the matrix 'A_big'
        for j=1:size(A_big,2)-L+1 %|     along the two directions rows/columns
            
            submatrix=A_big(i:i+L-1,j:j+L-1); % submatrix of 'A_big' of dimension LxL on which the window acts
            B(i,j)=sum(sum(submatrix.*g)); %| ==> the entry (i,j) of the post-filtered matrix is the sum of the entries
                                           %|     of the matrix corresponding to the entry per entry producta between
                                           %|     the submatrix of 'A_big' and the window
            
        end
    end
    
end

end
