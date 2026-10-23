function [modes,groups]=ct_cluster(z,centerHz,separationHz,labels)
%CT_CLUSTER Explicit grouping, or connected groups of adjacent centers.
% Automatic proximity clustering is an implementation choice; the paper
% does not specify a clustering algorithm. Harmonics require explicit labels.
k=size(z,1);
if nargin<4 || isempty(labels)
    labels=1:k;
    if separationHz>0 && k>0
        [s,ix]=sort(centerHz(:)); g=cumsum([1; diff(s)>separationHz]);
        labels(ix)=g;
    end
end
if numel(labels)~=k || any(~isfinite(labels))
    error('CT:Clusters','One finite cluster label per selected channel required.');
end
u=unique(labels,'stable'); modes=complex(zeros(numel(u),size(z,2)));
groups=cell(numel(u),1);
for j=1:numel(u)
    groups{j}=find(labels==u(j)); modes(j,:)=sum(z(groups{j},:),1);
end
end
