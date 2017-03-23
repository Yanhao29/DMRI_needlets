function separation_angle = separation_angle(pos1,pos2)
    pos1 = reshape(pos1,1,numel(pos1));
    pos2 = reshape(pos2,1,numel(pos2));
    separation_angle=acos(pos1*pos2');
end