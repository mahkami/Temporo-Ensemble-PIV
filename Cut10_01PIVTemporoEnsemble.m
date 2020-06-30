files = dir('*.tif');

load('frames.mat')
load('minimg.mat')

coeff_min = 1.1;

for i = 1: numel(files)

frames(:,:,i) = max(frames(:,:,i) - coeff_min*minimg,0);

end

height= size(frames(:,:,1),1);
width= size(frames(:,:,2),2);
searchy= 12;      % search range is 2*search+1
searchx= 6;
ox= 0;          % x-offset for search is 2*ox+1
oy= 0;          % y-offset for search is 2*oy+1
wxy= 2;         % spatial window size is 2*wxy+1
end_frame = 300;

coef_searchy = searchy/searchy;

%%

% loop over all grid points / pixels
nx= floor((width-2*(searchx+wxy+ox)+1)/(wxy+1));
ny= floor((height-2*(searchy+wxy+ox)+1)/(wxy+1));
xgrid= NaN(ny,nx);
ygrid= NaN(ny,nx);
ushift= NaN(ny,nx);
vshift= NaN(ny,nx);
peak= NaN(ny,nx);

ushift2= NaN(ny,nx);
vshift2= NaN(ny,nx);
peak2= NaN(ny,nx);

mm = 0;
% ii= 0;
% for i=1+search+wxy+ox:wxy+1:width-search-wxy-ox
%     ii= ii + 1;
%     disp (ii)
%     jj= 0;
%     for j=1+search+wxy+oy:wxy+1:height-search-wxy-oy
%         jj= jj + 1;
%         % first feature vector uses the same pixel in all frames
jj= 0;
for j=1+2*coef_searchy+wxy+oy:wxy+1:height-(2*searchy-2)-wxy-oy
    jj= jj + 1;
    disp (jj)
    ii= 0;
    for i=1+searchx+wxy+ox:wxy+1:width-searchx-wxy-ox
        ii= ii + 1;
        % first feature vector uses the same pixel in all frames
        av= frames(j-wxy-oy:j+wxy-oy,i-wxy-ox:i+wxy-ox,1:end_frame-1);
        %         av2 =uint8( mean(av,3));
        %         figure('name','av');surf(av2)
        for k=-searchx:searchx
            for l=-2*coef_searchy:(2*searchy-2)
                % 2nd feature vector uses the same pixel in all frames + 1
                bv= frames(j+l+oy-wxy:j+l+oy+wxy,i+k+ox-wxy:i+k+ox+wxy,2:end_frame);
                %                 bv2 = uint8(mean(bv,3));
                %                 figure();surf(bv2)
                % match criterion is correlation coefficent of feature vectors
                
                %                 gooddata = A(:,Col1)~=0 & A(:,Col2)~=0;
                %                 pearson = corr(A(gooddata,Col1),A(gooddata,Col2));
                
                
                %                 filt_indx = find(av(:)~=0 & bv(:)~=0);
                %                 filt_indx2 = find(av2(:)~=0 & bv2(:)~=0);
                %
                %                 corr(l+search+1,k+search+1)= corr2(av(filt_indx),bv(filt_indx));
                %                 corr_2(l+search+1,k+search+1)= corr2(av2(filt_indx2),bv2(filt_indx2));
                
                
                
                %                 ind_av = find(av(:)~=0);
                %                 ind_bv = find(bv(:)~=0);
                %                 ind_av2 = find(av2(:)~=0);
                %                 ind_bv2 = find(bv2(:)~=0);
                %
                %                 corr(l+search+1,k+search+1)= corr2(av(ind_av),bv(ind_bv));
                %                 corr_2(l+search+1,k+search+1)= corr2(av2(ind_av2),bv2(ind_bv2));
                
                
                corr(l+searchy+1,k+searchx+1)= corr2(av(:),bv(:));
                %                 corr_2(l+search+1,k+search+1)= corr2(av2(:),bv2(:));
                
                
            end
        end
        
        % search for shift with best match
        [maxval,maxind]= max(corr(:));
        [jmax,imax]= ind2sub (size(corr),maxind);
        
        %         [maxval2,maxind2]= max(corr_2(:));
        %         [jmax2,imax2]= ind2sub (size(corr_2),maxind2);
        
        
        
        % subpixel interpolation as usual
        [dx,dy,cmax]= GaussPeakFit (corr);
        
        %         [dx2,dy2,cmax2]= GaussPeakFit (corr_2);
        
        
        % Comparing two methods
        
        %         mm = mm+1;
        %         max_corr(mm,1) = maxval;
        %         max_corr(mm,2) = maxval2;
        %
        %         max_ind(mm,1) = maxind;
        %         max_ind(mm,2) = maxind2;
        %
        %         max_i(mm,1) = imax;
        %         max_i(mm,2) = imax2;
        %
        %         max_j(mm,1) = jmax;
        %         max_j(mm,2) = jmax2;
        
        % load into result arrays
        xgrid(jj,ii)= i;
        ygrid(jj,ii)= j;
        if cmax > 0.1 && cmax < 1 %imax>1 && jmax>1 && imax<=2*search && jmax<=2*search;
%         if maxval>0.5 && imax>1 && jmax>1 && imax<=2*search && jmax<=2*search

% if maxval>0.6 && maxval<1

               if isfinite(dy)
            ushift(jj,ii)= (imax-searchx-1) + dx; %+ 2*ox+1;
            vshift(jj,ii)= (jmax-searchy-1) + dy; %+ 2*oy+1;
            
            peak(jj,ii)= cmax;
               end
        else
            ushift(jj,ii)= NaN;
            vshift(jj,ii)= NaN;
            peak(jj,ii)= NaN;
        end
        
        
        %         if cmax2 > 0.1 && cmax2 < 1 %&& cmax2 > maxval2
        %             ushift2(jj,ii)= (imax2-search-1) + dx2; %+ 2*ox+1;
        %             vshift2(jj,ii)= (jmax2-search-1) + dy2; %+ 2*oy+1;
        %             peak2(jj,ii)= cmax2;
        %         else
        %             ushift2(jj,ii)= NaN;
        %             vshift2(jj,ii)= NaN;
        %             peak2(jj,ii)= NaN;
        %         end
        
        
%     figure (20)
%     quiver (xgrid,ygrid,ushift,vshift,5);
%     axis equal
%     drawnow
%         
%               
        
    end
    
%     figure (23)
%     quiver (xgrid,ygrid,ushift,vshift,5);
%     axis equal
%     drawnow
    
    
    %     figure (6)
    %     quiver (xgrid,ygrid,ushift2,vshift2,5);
    %     axis equal
    %     drawnow
    
    
    save InvWindows5_Shift25X13_10Cut_cmax01.mat xgrid ygrid ushift vshift peak searchx searchy ox oy wxy
    
end

save InvWindows5_Shift25X13_10Cut_cmax01.mat xgrid ygrid ushift vshift peak searchx searchy ox oy wxy


% 
% figure (23)
% quiver (xgrid,ygrid,ushift,vshift,5);
% axis equal
% drawnow
% 
% idx= find (isfinite(ushift(:)));
% figure (4)
% subplot (2,1,1)
% [counts,bins]= hist(ushift(idx),100);
% bar (bins,counts)
% xlabel ('x shifts [pxl]')
% ylabel ('counts')
% subplot (2,1,2)
% [counts,bins]= hist(vshift(idx),100);
% bar (bins,counts)
% xlabel ('y shifts [pxl]')
% ylabel ('counts')


% idx2= find (isfinite(ushift2(:)));
% figure (5)
% subplot (2,1,1)
% [counts,bins]= hist(ushift2(idx2),100);
% bar (bins,counts)
% xlabel ('x shifts [pxl]')
% ylabel ('counts')
% subplot (2,1,2)
% [counts,bins]= hist(vshift2(idx),100);
% bar (bins,counts)
% xlabel ('y shifts [pxl]')
% ylabel ('counts')

