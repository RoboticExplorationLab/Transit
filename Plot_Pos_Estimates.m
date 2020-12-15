%Loop over a range of velocity, frequency error estimates, and headings 
%and make a plot

%what ranges of inputs do we eval over?
minv=0;
maxv=10; %m/sec
min_freq=0.5;
max_freq=150;   %Hz
min_heading=0;
max_heading=pi/2;   %radians

%how many bins do we evaluate over
numvel_eval=10;
numfreq_eval=15;
numheading_eval=3;


Z=zeros(numfreq_eval,nx,numheading_eval);
[Vel,Freq,Head] = meshgrid(linspace(minv,maxv,nx),...
    logspace(log10(min_freq),log10(max_freq),numfreq_eval),...
    linspace(min_heading,max_heading,numheading_eval));
for i=1:nx
    for j=1:numfreq_eval
        for k=1:numheading_eval
            Z(j,i,k) = PositionEstimator(Vel(j,i,k),Head(j,i,k),Freq(j,i,k));
            %some random values are extremely bad; reject
            cnt=1;
            while Z(j,i,k)>2e6 && cnt<5
                Z(j,i,k) = PositionEstimator(Vel(j,i,k),Head(j,i,k),Freq(j,i,k));
            end
        end
    end
end

%spread of the area to be computed
r=floor(sqrt(numheading_eval));
c=ceil(numheading_eval/r);
for l=1:numheading_eval
    subplot(r,c,l);
    s = surf(Vel(:,:,l),Freq(:,:,l),Z(:,:,l)/1000);
    s.EdgeColor = 'none';
    colorbar
    xlabel('Unmodeled Tag Velocity (m/s)');
    ylabel('Tag Freq 1-sigma error (Hz)');
    title({'Position Estimate Error (km)',['Heading ',num2str(Head(1,1,l)),' rad']});
    set(gca, 'YScale', 'log');
    set(gca,'ColorScale','log')
    view(2); %view from the top
end
