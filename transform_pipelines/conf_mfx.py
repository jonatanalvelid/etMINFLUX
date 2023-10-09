def conf_mfx(coords_input, params, *args, **kwargs):
        "Pre-calibrated coordinate transform for etMINFLUX."

        conf_top_x_mon = params[0]
        conf_top_y_mon = params[1]
        conf_size_px_mon = params[2]
        conf_size_um = params[3]
        conf_size_px = params[4]

        conf_px_size_mon = conf_size_um/conf_size_px_mon

        x_um = coords_input[0]*conf_size_um/conf_size_px
        y_um = coords_input[1]*conf_size_um/conf_size_px
        x_um_corr = x_um - conf_size_um/2
        y_um_corr = y_um - conf_size_um/2

        x_px_mon = conf_top_x_mon + x_um/conf_px_size_mon
        y_px_mon = conf_top_y_mon + y_um/conf_px_size_mon

        return [x_px_mon, y_px_mon], [x_um_corr, y_um_corr], conf_px_size_mon