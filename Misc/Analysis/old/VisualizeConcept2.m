function VisualizeConcept2(concept_no,A,Z,Q,obj_names,fea_names)

concept_Z = Z(:,concept_no)==1;
concept_Q = Q(:,concept_no)==1;

figure
%Current concept's objects
subplot(2,2,2)
imagesc(A(concept_Z,:))
set(gca,'YTick',1:length(concept_Z));
set(gca,'YTickLabel',obj_names(concept_Z));
set(gca,'XTick',find(concept_Q));
set(gca,'XTickLabel',1:sum(concept_Q));


title(['Group: ' num2str(concept_no) ' - ' num2str(sum(concept_Z))  ' objects'])

%Current concept's features
subplot(2,2,4)
imagesc(A(:,concept_Q)')
set(gca,'YTick',1:length(concept_Q));
set(gca,'YTickLabel',fea_names(concept_Q));
set(gca,'XTick',find(concept_Z));
set(gca,'XTickLabel',obj_names(concept_Z));
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