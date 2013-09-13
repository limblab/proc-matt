%% Make plots of behavior metrics
html = strcat(html,['<div id="behavior"><h2>Behavior Metrics</h2>' ...
    '<table><tr> <td>epoch</td> <td>Time To Target</td> <td>Reaction Time</td> <td>Target Directions</td></tr>']);

for iEpoch = 1:length(epochs)
    html = strcat(html,['<tr>' ...
        '<td>' epochs{iEpoch} '</td> <td><img src="' figPath '\' epochs{iEpoch} '_behavior_time_to_target.png" width="' num2str(imgWidth+100) '"></td>' ...
        '<td><img src="' figPath '\' epochs{iEpoch} '_behavior_reaction_time.png" width="' num2str(imgWidth+100) '"></td>' ...
        '<td><img src="' figPath '\' epochs{iEpoch} '_behavior_target_direction.png" width="' num2str(imgWidth+100) '"></td>' ...
        '</tr>']);
end

html = strcat(html,'</table><br><a href="#header">back to top</a></div><hr>');