var layout = {
  title: "",
  height: Math.max(data.length * data[0]["num_cancers"] * 30, 500),
  xaxis: {
    side: "top",
    title: {
      text: `${miRNA} (RPM)`,
      standoff: 0
    },
    zeroline: false,
  },
  yaxis: {
  },
  margin: {
    r: 0
  },
  boxmode: "group"
}

var config = {responsive: true}

Plotly.newPlot("expression-plot", data, layout, config);
