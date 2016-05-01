function x = ransac(matchedPts1, matchedPts2)

% Repeat N times:
%   Draw s points uniformly at random
%   Fit line to these s points
%   Find inliers to this line among the remaining points 
%       (i.e., points whose distance from the line is less than t)
%   If there are d or more inliers, accept the line and refit using all inliers

if (length(matchedPts1)~=length(matchedPts2))
    error('Lists of matching points must be the same size');
end

p = 0.99;
e = 0.05;
s = ceil(length(matchedPts1)/2);
N = Inf;
t = 5;
% t = sqrt(3.84*std(std([matchedPts1,matchedPts2])));
x = zeros(size(matchedPts1,1),2);
N_inliers = zeros(size(matchedPts1,1),1);
i = 1;


while N>i
    try
        randPt = randi([1 length(matchedPts1)],s,1);
        displacementx = matchedPts1(randPt,1)-matchedPts2(randPt,1);
        displacementy = matchedPts1(randPt,2)-matchedPts2(randPt,2);
        b=[displacementx; displacementy];
        A1 = repmat([1 0],s,1);
        A2 = repmat([0 1],s,1);
        A = [reshape(A1',[],1),reshape(A2',[],1)];
        x(i,:) = A \ b;
        remainder = setdiff(1:length(matchedPts1), randPt);
        translation1 = matchedPts1(remainder,:) - repmat(x(i,:),length(remainder),1);
        inliers = find(pdist2(translation1,matchedPts2(remainder,:)) < t);
        N_inliers(i) = length(inliers);

        % Adaptively determine number of samples
        e = 1 - length(inliers)/size(matchedPts1,1);
        N = log(1-p)/log(1-(1-e).^s);
        i = i + 1;
    catch
        x = zeros(size(matchedPts1,1),2);
        break;
    end
end

[~,idx] = max(N_inliers);
x = x(idx);