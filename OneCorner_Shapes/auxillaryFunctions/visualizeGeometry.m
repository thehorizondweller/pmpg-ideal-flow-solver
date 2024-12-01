function visualizeGeometry(z)
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    markerSize = 4; % Size of markers
    markerColor = 'b'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(z),imag(z),markerSize,markerColor,'filled',markerShape);
    xlabel("X");
    ylabel("Y");
    title("Physical Z Domain");
    axis equal;
    hold off;
    outputFile = fullfile('./Results','01_UserGeometry.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 1: Geometry Generation Complete.");
end