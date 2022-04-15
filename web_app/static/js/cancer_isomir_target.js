var layout = {
  title: `Spearman's r = ${spearman_corr}, <i>p</i>-value = ${spearman_p_value}`,
  height: 600,
  yaxis: {
    title: {
      text: `${target} (log<sub>2</sub> TPM)`
    }
  },
  xaxis: {
    title: {
      text: `${isomir} (log<sub>2</sub> RPM)`
    }
  }
}

var config = {responsive: true}
data[0]["showlegend"] = false;
data[0]["marker"] = {"size": 8};
//data[0]["opacity"] = 0.999;
data[1]["showlegend"] = false;

Plotly.newPlot(`expression-plot-${cancer}`, data, layout, config);
