data[0]["box"] = {visible: false}
data[0]["side"] = "positive"
data[0]["meanline"] = {visible: true}
data[0]["points"] = false
data[0]["spanmode"] = "hard"

var layout = {
  title: "",
  yaxis: {
    showticklabels: false
  },
  xaxis: {
    title: {
      text: `${molecule} (log<sub>2</sub> RPM)`,
    },
    zeroline: false,
  }
}

var config = {responsive: true}

Plotly.newPlot("expression-plot", data, layout, config);
