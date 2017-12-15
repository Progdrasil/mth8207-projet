function saveImg(fig, title)
	fig_path = fullfile('graphs', title);
	saveas(fig, fig_path, 'png');
end
