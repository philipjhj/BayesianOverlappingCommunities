function VisualizeConcept(concept_no,A,Z,Q,obj_names,fea_names)

concept_Z = Z(:,concept_no)==1;
concept_Q = Q(:,concept_no)==1;

figure
%Current concept's objects
subplot(2,2,2)
[sorted_obj_concepts, idx1_obj, idx1_fea] = sort_matrix(A(concept_Z,:),1);
imagesc(sorted_obj_concepts)
set(gca,'YTick',1:length(concept_Z));
names_for_obj = obj_names(concept_Z);
set(gca,'YTickLabel',names_for_obj(idx1_obj));
title(['Group: ' num2str(concept_no) ' - ' num2str(sum(concept_Z))  ' objects'])

%Current concept's features
subplot(2,2,4)
[sorted_fea_concepts, idx2_obj, idx2_fea] = sort_matrix(A(:,concept_Q),1);
imagesc(sorted_fea_concepts')
set(gca,'YTick',1:length(concept_Q));
names_for_fea = fea_names(concept_Q);
set(gca,'YTickLabel',names_for_fea(idx2_fea));
title(['Group: ' num2str(concept_no) ' - ' num2str(sum(concept_Q)) ' features'])

%Combination of current concept's objects and features
subplot(2,2,[1 3])
Concept_A = A;
Concept_A(concept_Z,concept_Q) = A(concept_Z,concept_Q)*2;

imagesc(Concept_A')
set(gca,'YTick',find(concept_Q));
set(gca,'YTickLabel',fea_names(concept_Q));
set(gca,'XTick',find(concept_Z));
set(gca,'XTickLabel',obj_names(concept_Z));
title(['Group: ' num2str(concept_no) ' - ' num2str(sum(concept_Z)) ' objects - ' num2str(sum(concept_Q)) ' features'])

end