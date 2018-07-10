// Author Max Briel
// Known Gravitational wave events
var events_list = [
    "GW150914",
    "LVT151012",
    "GW151226",
    "GW170104",
    "GW170608",
    "GW170814",
    "GW170817"
];


// Set the width of the page
var contdiv = document.getElementById("container");
var width = contdiv.clientWidth;
var height = contdiv.clientHeight;

// Create an svg in the images div
var svg = d3.select("#images").select("svg");

svg.attr("viewBox","5,1.5,100,50");
//
// attr("width", "100%")
//     .attr("height", "60%");

// add the event image plots to the svg
svg.selectAll("image").data(events_list).enter()
    .append("image")
        .attr("id", function(d) {return "image_event_"+d;})
        .attr("width", "100%")
        .attr("height","100%")
        .style("opacity", 0)
        .attr("xlink:href", function(d) {return "events/"+d+"/"+d+".png";});

// Add the background graticule to the svg on top
svg.insert("image", ":first-child")
        .attr("id", "background")
        .attr("width", "100%")
        .attr("height","100%")
        .attr("xlink:href", "events/graticule.png");

// selects all divs in the event_selector div
var columns = d3.select("#images").select("#event_selector").selectAll("div");
// add 2 divs for events per column
columns.each(function(d,i) {
    var col = d3.select(this);
    col.append('div')
        .attr('class',"event");
    col.append("div")
        .attr('class', 'event');
});

// add a checkbox, an label, and listener for each event div
var events = d3.select("#event_selector").selectAll(".event");
events.each(function(d, i){
    var event = d3.select(this);
    event.attr("id","event_"+events_list[i]);
    event.append("label")
        .append('input')
            .attr("type", "checkbox")
            .attr("id", "checkbox_"+events_list[i]);
    event.select("label")
        .append("text")
        .text(events_list[i]);
    event.on("change", update);
    if (i > events_list.length-1){
        event.remove();
    }
});

// update function to change if the image and data of the event is shown
function update() {
    var event = d3.select(this);
    // hide the image of the event
    var event_image = d3.select("#image_"+event.attr("id"));
    var opacity = event_image.style("opacity");
    if (opacity == 0){
        event_image.style("opacity", 1);
    }
    else {
        event_image.style("opacity", 0);
    }
    // hide the data of the event
    var event_data = d3.select("#data_"+event.attr("id"));
    var opacity = event_data.style("display");
    if (opacity == "inline"){
        event_data.style("display", "none");
    }
    else {
        event_data.style("display", "inline");
    }
};


d3.select("#data_showcase").selectAll("div").data(event_data).enter()
    .append("div")
        .attr("class","data_container")
        .attr("id",function (d) {return "data_event_"+d.name})
        .style("display", "none");

d3.select("#data_showcase").selectAll(".data_container").each(
function (d,i) {
    var data_container = d3.select(this);

    // Add a table containing event information to the conatiner
    data_container.append("div")
        .attr("class", "title")
        .text(d.name)
        .style("color","#171717")
        .style("font-weight", "bold")
        .style("text-align","center");
    var data = data_container.append("div").attr("id", "data_"+d.name)
    var left_column = data.append("div").attr("class", "left_column");

    var table = left_column.append("table");
    var thead = table.append("thead");
    var tbody = table.append("tbody");
    // Create header of table
    thead.append("tr")
        .selectAll('th')
        .data(["Quantity", "Value"]).enter()
        .append("th")
            .text(function(header, i) {
                if (i ==1){
                    d3.select(this).style("text-align", "right");
                }; return header;});

    // Add the rows with data
    var rows = tbody.selectAll("tr")
        .data(function(j) {
            var keys = Object.keys(j);
            var values = Object.values(j);
            out = [];
            for (i=1; i<keys.length;i++ ){
                out.push([keys[i], values[i]]);
            }
            return out;
        }).enter()
            .append("tr");

        rows.selectAll("td").data(function(d) {return d;}).enter()
            .append("td").text(function(d, i) {
                if (i==1){
                    d3.select(this).style("text-align","right");}; return d;});

    // Add spectogram
    var sound = left_column.append("div").attr("width", "100%");
    sound.append("img")
        .attr("src", "events/"+d.name+"/spectrogram.png")
        .attr("width", "95%")
        .style("margin","15px 0px 0px 0px");

    // Add audio gile to play
    sound.append("p")
        .text("Click play to hear the audio:")
        .style("text-align", "center")
        .style("margin", "20px 10px 5px 10px");
    sound.append("audio")
        .attr("id", "player"+d.name)
        .attr("src","events/"+d.name+"/sound.wav");

    // Add the play button and connect it to the audio
    sound.append("div")
        .style("text-align", "center")
        .append("button")
            .on("click",function() {
                document.getElementById("player"+d.name).play();
            })
            .text("Play");

    // Add plots of the event
    var data_images = data.append("div").attr("class","data_images");
    data_images.append("img")
        .attr("src", "events/"+d.name+"/waveform.png")
        .attr("width","100%");
    data_images.append("img")
        .attr("src", "events/"+d.name+"/ASD.png")
        .attr("width","100%");
});


// Set one GW event to visible and check its checkbox
var event = d3.select(".event");
event.select("input").property("checked", true);
d3.select("#image_"+event.attr("id")).style("opacity", 1);
d3.select("#data_"+event.attr("id")).style("display", "inline");

// Add informational text
var text = d3.select("#images").append("div");
text.attr("id", "exp_text").html(
    "<h2> Gravitational Wave Detections Summarised</h2> "
    + "<p> The Mollweide projection map show the probabilities over the whole "
    + "sky. The most probable location and region are drawn. During the "
    + "first runs no distinction could be made between a position above or "
    + "below the detectors. Therefore, mulitple disconnected regions are "
    + "present. </p> On the right the waveform and amplitude spectral density "
    + "(log-log scale) from Hanford 1 during the event is shown. "
    + "Furthermore, a spectrogram of the event gives insight on the "
    + "evolution of the frequency. This can be heard in the audio file."
    + "<p>Detected Gravitational wave events up to May 20th 2018 "
    + "with data from the LIGO H1 detector.<br/>"
    + "An mBriel Project. Original data for the events can be found at the "
    + "<a href='https://losc.ligo.org/events/'>LIGO website</a>. <br/>"
).style("margin-top", "-20px");
