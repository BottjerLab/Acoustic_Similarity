% type the syllables
num2cell(clusterIdxs(:,end)); [typedDRsylls.type] = ans{:};

N = numel(typedDRsylls);
conformist = zeros(2,N); % zero is most conformist, 
takenIdxs = [typedDRsylls.type];
sqcosim = squareform(distMats.cosim);
for ii = 1:N
    thisType = takenIdxs(ii);
    isAlsoType = (takenIdxs == thisType);
    isNotType = ~isAlsoType;
    isAlsoType(ii) = false;
    conformist(1,ii) = mean(sqcosim(ii,isAlsoType))/mean(sqcosim(ii,isAlsoType | isNotType));
    conformist(2,ii) = var(sqcosim(ii,isAlsoType))/var(sqcosim(ii,isAlsoType | isNotType));
end
figure(10)
subplot(2,1,1)
boxplot(conformist(1,:)',takenIdxs')
subplot(2,1,2)
boxplot(conformist(2,:)',takenIdxs')

figure(11)
nTypes = max(takenIdxs);
cols = jet(nTypes);
for ii = 1:nTypes
    hh(ii) = plot(conformist(1,takenIdxs == ii), conformist(2,takenIdxs == ii), '.', 'Color', cols(ii,:));
    hold on;
end

legLabels = arrayfun(@num2str,1:nTypes,'UniformOutput',false)';
legend(legLabels);
hold off;