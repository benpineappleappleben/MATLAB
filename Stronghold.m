function c = Stronghold(Theta1,x1,z1,Theta2,x2,z2)
% Stronghold() has six inputs (in order): Theta1, x1, z1, Theta2, x2, z2.
% x1 and z1 are the player's coordinates X and Z at the first Eye of Ender 
% throw, and Theta1 is the Minecraft displayed angle on the X-Z plane for 
% that same throw. Similarly, x2 and z2 are the player's coordinates X and 
% Z at the second Eye of Ender throw, and Theta2 is the Minecraft displayed
% angle on the X-Z plane for the same throw. All six of these variables may
% be found in the debugging menu, accessible by pressing F3 in-game. Here
% is what the Minecraft debugging menu looks like and where to find each of
% the six variables:
%
% Minecraft 1.17 (Version)
% ... fps t: ......
% ... @ ... ms ticks, ... tx, ... rx
% C: ...
% E: ...
% P: ...
% Chunks[C] W:...
% Chunks[S] W:...
% minecraft:overworld ...
% 
% XYZ: <strong>x1 or x2</strong> / ... / <strong>z1 or z2</strong>
% Block: ...
% Chunk: ...
% Facing: ... (Towards ...) (<strong>Theta1 or Theta 2</strong> / ...)
% <strong>You don't need anything after the above line </strong>
%
% Stronghold() uses these six parameters to triangulate the coordinates of
% the nearest stronghold.
%
% For the most accurate coordinates of the stronghold, make sure that there
% is plenty of distance between both Eye of Ender throws. Around 200 blocks
% is sufficient, although the function will approximate the coordinates for
% the stronghold regardless of how much distance there is between the two
% throws (distance cannot equal zero). In addition, be sure to point the
% crosshairs directly at the center of the Eye of Ender when recording the
% angle.
%
% The coordinates X and Z of the stronghold are displayed as the matrix
% [x,z]. Make sure that you type "[x,z] = Stronghold(...)" as opposed to
% "Stronghold(...)". The latter will only return the x-coordinate of the
% stronghold. In order to get both coordinates, you will have to enter in
% the former syntax.
% The syntax for the Stronghold() function itself is:
%
% [x,z] = Stronghold(Theta1,x1,z1,Theta2,x2,z2)
Theta1rad = Theta1*(pi/180);
Theta2rad = Theta2*(pi/180);
x = ((x1-x2+z1*tan(Theta1rad)-z2*tan(Theta2rad))/(tan(Theta1rad)-tan(Theta2rad)))*tan(Theta2rad) - x2 - z2*tan(Theta2rad);
x = -1 * x;
z = ((x1-x2+z1*tan(Theta1rad)-z2*tan(Theta2rad))/(tan(Theta1rad)-tan(Theta2rad)));
c(1) = x;
c(2) = z;
end