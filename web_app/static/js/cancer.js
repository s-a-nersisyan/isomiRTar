var nodes = new vis.DataSet(nodes);
var edges = new vis.DataSet(edges);
// create a network
var container = document.getElementById("network");

// provide the data in the vis format
var data = {
  nodes: nodes,
  edges: edges
};
var options = {
  width: "100%",
  height: "600px",
  nodes: {
    shape: "dot",
  },
  edges: {
    arrows: {
      to: {
        enabled: true,
        type: "bar"
      }
    },
    smooth: true
  },
  layout: {
    randomSeed: 2,
    improvedLayout: false,
  },
  physics: {
    enabled: true,
    solver: "forceAtlas2Based",
    minVelocity: 1,
    stabilization: {
      enabled: true,
      iterations: 1000,
      updateInterval: 10
    }
  },
  interaction: { multiselect: true}
};

// initialize your network!
var network = new vis.Network(container, data, options);

network.on("stabilizationProgress", function (params) {
  var perc = Math.floor(params.iterations / params.total * 100);
  document.getElementById("loader").innerHTML = `Loading the network<br>${perc}%`;
});
network.once("stabilizationIterationsDone", function () {
  // document.getElementById("loader").style.opacity = 0;
  // really clean the dom element
  document.getElementById("loader").style.display = "none";
  network.setOptions({physics: false});
});
