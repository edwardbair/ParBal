function S=extractInd(S,ind)
%extract inputs for time of day when Tsfc image is from
%input - S struct to extract from, ind - index for struct
fn=fieldnames(S);
    for j=1:length(fn)
        if ~strncmp(fn{j},'date',4)
            S.(fn{j})=S.(fn{j})(:,:,ind);
        else
            S.(fn{j})=S.(fn{j})(ind);
        end
    end