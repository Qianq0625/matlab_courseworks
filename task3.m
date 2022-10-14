%% Initialise the camera
% 
% clear; clc;
% K1 = [[1.02408777e+03,0.,6.01806702e+02];[0.,1.02389478e+03,5.08131683e+02];[0.,0.,1.]];
% D1 = [-2.52257544e-03,4.38565714e-03,0.,0.,9.21500759e-05];
% K2 = [[1.02419836e+03,0.,6.96750427e+02];[0.,1.02398749e+03,5.07494263e+02];[0.,0.,1.]];
% D2 = [-3.26306466e-03,5.70008671e-03,0.,0.,7.57322850e-05];
% RR = [[9.99999404e-01,-1.24090354e-06,-1.11142127e-03];[1.29263049e-06,1.,4.65405355e-05];[1.11142127e-03,-4.65419434e-05,9.99999404e-01]];
% T = [[-4.34818459e+00];[2.83603016e-02];[-9.00963729e-04]];
% 
% image_left = imread('.\Stereo Data\Stereo Data\Left_Image.png');
% image_right = imread('.\Stereo Data\Stereo Data\Right_Image.png');
% 
% camera_left = cameraParameters('IntrinsicMatrix',K1','RadialDistortion',[D1(1),D1(2),D1(5)],'TangentialDistortion',[D1(3),D1(4)]);
% camera_right = cameraParameters('IntrinsicMatrix',K2','RadialDistortion',[D2(1),D2(2),D2(5)],'TangentialDistortion',[D2(3),D2(4)]);

image_undistorted_left = undistortImage(image_left,camera_left);
image_undistorted_right = undistortImage(image_right,camera_right);
IntrinsicMatrix = [1.02408777e+03 0 6.01806702e+02; 0 1.02389478e+03 5.08131683e+02; 0 0 1];
IntrinsicMatrix = IntrinsicMatrix';
radialDistortion = [-2.52257544e-03 4.38565714e-03 9.21500759e-05]; 
TangentialDistortion=[0 0];
camera_left = cameraParameters('IntrinsicMatrix',IntrinsicMatrix,'RadialDistortion',radialDistortion,'TangentialDistortion',TangentialDistortion); 
IntrinsicMatrix = [1.02419836e+03 0 6.96750427e+02; 0 1.02398749e+03 5.07494263e+02; 0 0 1];
IntrinsicMatrix = IntrinsicMatrix';
radialDistortion = [-3.26306466e-03 5.70008671e-03 7.57322850e-05]; 
TangentialDistortion=[0 0];

camera_right = cameraParameters('IntrinsicMatrix',IntrinsicMatrix,'RadialDistortion',radialDistortion,'TangentialDistortion',TangentialDistortion); 

%% a. Detect salient features and plot the detected features on the provided pair of frames.

image_gray_left = im2gray(image_left);
image_gray_right = im2gray(image_right);

image_undistorted_gray_left = im2gray(image_undistorted_left);
image_undistorted_gray_right = im2gray(image_undistorted_right);

points_left = detectSURFFeatures(image_gray_left, 'MetricThreshold', 2000);
points_right = detectSURFFeatures(image_gray_right, 'MetricThreshold', 2000);

points_undistorted_left = detectSURFFeatures(image_undistorted_gray_left, 'MetricThreshold', 2000);
points_undistorted_right = detectSURFFeatures(image_undistorted_gray_right, 'MetricThreshold', 2000);


% figure;
% plot(points_left);
% title('Points Left');
% plot(points_right);
% title('Points Right');
% plot(points_undistorted_left);
% title('Points Undistorted Left');
% plot(points_undistorted_right);
% title('Points Undistorted Right');

image_points_left = imread("output_images\left.png");
image_points_right = imread("output_images\right.png");
figure;imshowpair(image_points_left,image_points_right,'montage');
title('100 strongest SURF features of Left frame         100 strongest SURF features of Right frame');

figure;imshow(image_left);hold on;
plot(selectStrongest(points_left, 100));
title('100 strongest SURF features in the left frame');

figure;imshow(image_right);hold on;
plot(selectStrongest(points_right, 100));
title('100 strongest SURF features in the right frame');
%% b. Find corresponding features between the two frames. Create a composite image from the two frames and demonstrate the matched points.

[image_features_left, points_valid_left] = extractFeatures(image_gray_left, points_left);
[image_features_right, points_valid_right] = extractFeatures(image_gray_right, points_right);

[image_undistorted_features_left, points_undistorted_valid_left] = extractFeatures(image_undistorted_gray_left, points_undistorted_left);
[image_undistorted_features_right, points_undistorted_valid_right] = extractFeatures(image_undistorted_gray_right, points_undistorted_right);

points_pairs_index = matchFeatures(image_features_left, image_features_right, 'Metric', 'SAD', ...
  'MatchThreshold', 5);

points_undistorted_pairs_index = matchFeatures(image_undistorted_features_left, image_undistorted_features_right, 'Metric', 'SAD', ...
  'MatchThreshold', 5);

points_matched_left = points_valid_left(points_pairs_index(:,1),:);
points_matched_right = points_valid_right(points_pairs_index(:,2),:);

points_undistorted_matched_left = points_undistorted_valid_left(points_undistorted_pairs_index(:,1),:);
points_undistorted_matched_right = points_undistorted_valid_right(points_undistorted_pairs_index(:,2),:);

figure;
showMatchedFeatures(image_left, image_right, points_matched_left, points_matched_right);
legend('Putatively matched points in the left frame(70)', 'Putatively matched points in the right frame(70)');
figure;
showMatchedFeatures(image_left, image_right, points_undistorted_matched_left, points_undistorted_matched_right);
legend('Putatively matched undistorted points in the left frame(70)', 'Putatively matched undistorted points in the right frame(70)');
%% c. Use the matched features to estimate the fundamental matrix between the two images. 

[fMatrix, epipolarInliers, status] = estimateFundamentalMatrix(...
  points_matched_left, points_matched_right, 'Method', 'RANSAC', ...
  'NumTrials', 10000, 'DistanceThreshold', 0.1, 'Confidence', 99.99);

[fMatrix_undistorted, epipolarInliers_undistorted, status_undistorted] = estimateFundamentalMatrix(...
  points_undistorted_matched_left, points_undistorted_matched_right, 'Method', 'RANSAC', ...
  'NumTrials', 10000, 'DistanceThreshold', 0.1, 'Confidence', 99.99);

if status ~= 0 || isEpipoleInImage(fMatrix, size(image_left)) ...
  || isEpipoleInImage(fMatrix', size(image_right))
  error(['Either not enough matching points were found or '...
         'the epipoles are inside the images. You may need to '...
         'inspect and improve the quality of detected features ',...
         'and/or improve the quality of your images.']);
end

if status_undistorted ~= 0 || isEpipoleInImage(fMatrix_undistorted, size(image_undistorted_left)) ...
  || isEpipoleInImage(fMatrix_undistorted', size(image_undistorted_right))
  error(['Either not enough matching points were found or '...
         'the epipoles are inside the images. You may need to '...
         'inspect and improve the quality of detected features ',...
         'and/or improve the quality of your images.']);
end

% F = e’*（P'P+)
F = inv(camera_right.IntrinsicMatrix)*skew(T)*RR*inv(camera_left.IntrinsicMatrix');
%solve
value1 = [(points_matched_right(1).Location),1]*F*[(points_matched_left(1).Location),1]';
% value 是两种方法误差值，基于x'Fx这个计算的，参照ppt topic 5

value2 = [(points_matched_right(1).Location),1]*fMatrix*[(points_matched_left(1).Location),1]';
value_undistorted2 = [(points_undistorted_matched_right(1).Location),1]*fMatrix_undistorted*[(points_undistorted_matched_left(1).Location),1]';
%% d. Find the correctly matched points that meet the epipolar constraint and illustrate these matches. Briefly explain how these matches have been identified.

[E, Inliers] = estimateFundamentalMatrix(points_undistorted_matched_left, points_undistorted_matched_right,'Method', 'RANSAC', ...
  'NumTrials', 10000, 'DistanceThreshold', 0.1, 'Confidence', 99.99);

inlierPoints1 = points_undistorted_matched_left(Inliers,:);
inlierPoints2 = points_undistorted_matched_right(Inliers,:);
figure
showMatchedFeatures(image_undistorted_left,image_undistorted_right,inlierPoints1,inlierPoints2);
title('Inlier Matches')

[orient, loc] = relativeCameraPose(E, camera_left, inlierPoints1, inlierPoints2);


% inlierPoints1 = points_matched_left(epipolarInliers, :);
% inlierPoints2 = points_matched_right(epipolarInliers, :);
% 
% figure;
% showMatchedFeatures(image_left, image_right, inlierPoints1, inlierPoints2);
% legend('Inlier points in the left frame', 'Inlier points in the right frame');
% 
% t_x=skew(T);
% E=t_x*RR;
% 
% F=((K2.')^-1)*E*((K1)^-1);
% 
% MP1=points_matched_left.Location;
% MP2=points_matched_right.Location;

%%
% Detect dense feature points. Use an ROI to exclude points close to the
% image edges.
[orient, loc] = relativeCameraPose(E, camera_left, camera_right, inlierPoints1, inlierPoints2);
% orient = orient(:,:,1);
% loc = loc(1,:);
% Detect dense feature points. Use an ROI to exclude points close to the
% image edges.
roi = [30, 30, size(image_undistorted_left, 2) - 30, size(image_undistorted_left, 1) - 30];
imagePoints1 = detectMinEigenFeatures(im2gray(image_undistorted_left), 'ROI', roi, ...
    'MinQuality', 0.0001);

% Create the point tracker
tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

% Initialize the point tracker
imagePoints1 = imagePoints1.Location;
initialize(tracker, imagePoints1, image_left);

% Track the points
[imagePoints2, validIdx] = step(tracker, image_undistorted_right);
matchedPoints1 = imagePoints1(validIdx, :);
matchedPoints2 = imagePoints2(validIdx, :);
figure
showMatchedFeatures(image_undistorted_left, image_undistorted_right, matchedPoints1, matchedPoints2);
title('MinEigen 0.0001 Features');

% Compute the camera matrices for each position of the camera
% The first camera is at the origin looking along the Z-axis. Thus, its
% rotation matrix is identity, and its translation vector is 0.
camMatrix1 = cameraMatrix(image_undistorted_left, eye(3), [0 0 0]);

% Compute extrinsics of the second camera
[R, t] = cameraPoseToExtrinsics(orient, loc);
% R = R_vec;
% t = T_vec';
camMatrix2 = cameraMatrix(image_undistorted_right, R, t);
% [orient, loc] = extrinsicsToCameraPose(R,t);
% Compute the 3-D points
points3D = triangulate(matchedPoints1, matchedPoints2, camMatrix1, camMatrix2);

% Get the color of each reconstructed point
numPixels = size(image_undistorted_left, 1) * size(image_undistorted_left, 2);
allColors = reshape(image_undistorted_left, [numPixels, 3]);
colorIdx = sub2ind([size(image_undistorted_left, 1), size(image_undistorted_left, 2)], round(matchedPoints1(:,2)), ...
    round(matchedPoints1(:, 1)));
color = allColors(colorIdx, :);

% Create the point cloud
ptCloud = pointCloud(points3D, 'Color', color);
% Visualize the camera locations and orientations
cameraSize = 0.5;
figure
% plotCamera('Size', cameraSize, 'Color', 'r', 'Label', '1', 'Opacity', 0);
hold on
grid on
% plotCamera('Location', loc, 'Orientation', orient, 'Size', cameraSize, ...
%     'Color', 'b', 'Label', '2', 'Opacity', 0);

% Visualize the point cloud
pcshow(ptCloud, 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
    'MarkerSize', 45);

% Rotate and zoom the plot
camorbit(0, -30);
camzoom(1.5);

% Label the axes
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis')

title('Up to Scale Reconstruction of the Scene');

% Detect the globe
globe = pcfitcylinder(ptCloud, 0.01);

% Display the surface of the globe
plot(globe);
title('Estimated Location and Size of the Globe');
hold off