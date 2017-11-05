function [ x_pos, y_pos, values ] = find2dpeaks( spectre )

[x y] = size(spectre);
x_pos = [];
y_pos = [];
h = zeros(x,y);
l = zeros(x,y);
for k=1:x
    [~, peaks] = findpeaks(spectre(k,:));
    h(k,peaks) = 1;
end
for k=1:y
    [~, peaks] = findpeaks(spectre(:,k));
    l(peaks,k) = 1;
end

[x_pos, y_pos] = find((h==1 & l==1));

values = diag(spectre(x_pos,y_pos));
[values, index] = sort(values,1,'descend');
x_pos = x_pos(index);
y_pos = y_pos(index);

end